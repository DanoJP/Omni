#!/usr/bin/env python3
import RPi.GPIO as GPIO
from adafruit_pca9685 import PCA9685
import busio
from board import SCL, SDA
import rospy
from geometry_msgs.msg import Twist
import numpy as np
from kinematics import E  # Importar matriz E

# Configuración de los pines GPIO
GPIO.setmode(GPIO.BCM)
motor_pins = {
    "M1": (23, 22),
    "M2": (27, 18),
    "M3": (17, 4)
}
for pins in motor_pins.values():
    GPIO.setup(pins, GPIO.OUT)
    GPIO.output(pins, False)

# Configuración del PCA9685
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1000  # Frecuencia PWM (1 kHz)
MAX_VOLTAGE = 12  # Voltaje máximo del motor

def set_motor_speed(motor_pins, speed, channel):
    """Controla un motor con pines (IN1, IN2) y velocidad [-MAX_VOLTAGE, MAX_VOLTAGE]."""
    if speed >= 0:
        GPIO.output(motor_pins[0], False)
        GPIO.output(motor_pins[1], True)
    else:
        GPIO.output(motor_pins[0], True)
        GPIO.output(motor_pins[1], False)

    # Convertir velocidad a duty cycle PWM (0 - 65535)
    duty_cycle = int(min(abs(speed) * 65535 / MAX_VOLTAGE, 65535))
    pca.channels[channel].duty_cycle = duty_cycle

def cmd_vel_callback(msg):
    """Convierte Twist en velocidades de motores usando cinemática inversa."""
    vx, vy, omega = msg.linear.x, msg.linear.y, msg.angular.z
    motor_speeds = np.dot(E, [vx, vy, omega])

    # Enviar velocidades a los motores
    set_motor_speed(motor_pins["M1"], motor_speeds[0], 0)
    set_motor_speed(motor_pins["M2"], motor_speeds[1], 1)
    set_motor_speed(motor_pins["M3"], motor_speeds[2], 2)

if __name__ == '__main__':
    rospy.init_node('motor_control')
    rospy.Subscriber('/cmd_vel', Twist, cmd_vel_callback)
    rospy.spin()
    GPIO.cleanup()  # Limpiar GPIO al cerrar
