#!/usr/bin/env python3
import rospy
import numpy as np
from geometry_msgs.msg import Twist, PoseStamped
from tf.transformations import euler_from_quaternion
from kinematics import E  # Importar la matriz E

# Ganancias del controlador PD
Kp = np.diag([2.0, 2.0, 1.5])  # Proporcional (x, y, theta)
Kd = np.diag([0.5, 0.5, 0.3])  # Derivativo (x, y, theta)

class RobotController:
    def __init__(self):
        self.current_pose = np.zeros(3)  # [x, y, theta]
        self.last_error = np.zeros(3)
        self.last_time = rospy.get_time()

        # Suscriptores y publicadores
        rospy.Subscriber('/robot_position', PoseStamped, self.pose_callback)
        rospy.Subscriber('/setpoint', Twist, self.setpoint_callback)
        self.cmd_vel_pub = rospy.Publisher('/cmd_vel', Twist, queue_size=10)

        # Variables de la trayectoria deseada
        self.desired_pose = np.zeros(3)  # [x_d, y_d, theta_d]
        self.desired_vel = np.zeros(3)   # [vx_d, vy_d, omega_d]

    def pose_callback(self, msg):
        """Actualiza la posición actual del robot desde OptiTrack."""
        self.current_pose[0] = msg.pose.position.x
        self.current_pose[1] = msg.pose.position.y
        orientation_q = msg.pose.orientation
        _, _, self.current_pose[2] = euler_from_quaternion([orientation_q.x, orientation_q.y, orientation_q.z, orientation_q.w])

    def setpoint_callback(self, msg):
        """Actualiza el setpoint de la trayectoria deseada."""
        self.desired_pose[0] = msg.linear.x
        self.desired_pose[1] = msg.linear.y
        self.desired_vel[0] = msg.angular.x  # vx_d
        self.desired_vel[1] = msg.angular.y  # vy_d
        self.desired_vel[2] = msg.angular.z  # omega_d

    def compute_pd_control(self):
        """Calcula el comando de velocidad usando control PD."""
        current_time = rospy.get_time()
        dt = max(current_time - self.last_time, 0.01)  # Evitar dt = 0

        # Error de posición y velocidad
        error_pos = self.desired_pose - self.current_pose
        error_vel = (error_pos - self.last_error) / dt

        # Control PD
        control = np.dot(Kp, error_pos) + np.dot(Kd, error_vel) + self.desired_vel

        # Guardar error y tiempo para la siguiente iteración
        self.last_error = error_pos
        self.last_time = current_time

        return control

    def run(self):
        rate = rospy.Rate(20)  # 20 Hz
        while not rospy.is_shutdown():
            control = self.compute_pd_control()

            # Publicar comando de velocidad
            cmd_vel_msg = Twist()
            cmd_vel_msg.linear.x = control[0]
            cmd_vel_msg.linear.y = control[1]
            cmd_vel_msg.angular.z = control[2]
            self.cmd_vel_pub.publish(cmd_vel_msg)
            
            rate.sleep()

if __name__ == '__main__':
    rospy.init_node('robot_control')
    controller = RobotController()
    controller.run()
