#!/usr/bin/env python3
import rospy
import numpy as np
from geometry_msgs.msg import Twist

class TrajectoryGenerator:
    def __init__(self):
        self.pub = rospy.Publisher('/setpoint', Twist, queue_size=10)
        self.radius = 0.5  # Radio de la trayectoria (m)
        self.omega = 0.5   # Velocidad angular (rad/s)
        self.rate = rospy.Rate(10)  # 10 Hz

    def generate_circle(self):
        t = 0.0
        while not rospy.is_shutdown():
            # Posición deseada
            x_d = self.radius * np.cos(self.omega * t)
            y_d = self.radius * np.sin(self.omega * t)
            
            # Velocidad deseada (derivada de la posición)
            vx_d = -self.radius * self.omega * np.sin(self.omega * t)
            vy_d = self.radius * self.omega * np.cos(self.omega * t)
            
            # Crear mensaje
            setpoint = Twist()
            setpoint.linear.x = x_d
            setpoint.linear.y = y_d
            setpoint.angular.x = vx_d  # Velocidad deseada en X
            setpoint.angular.y = vy_d  # Velocidad deseada en Y
            setpoint.angular.z = 0.0   # Velocidad angular deseada
            
            self.pub.publish(setpoint)
            t += 0.1  # Incrementar tiempo (dt = 1/rate)
            self.rate.sleep()

if __name__ == '__main__':
    rospy.init_node('trajectory_generator')
    generator = TrajectoryGenerator()
    generator.generate_circle()