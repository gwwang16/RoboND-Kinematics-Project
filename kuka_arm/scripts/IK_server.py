#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, lambdify
from sympy.matrices import Matrix
import pickle, os
import time

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []

        def DH_transform(q, d, alpha, a): 
            '''The modified DH convention transform matrix
            alpha: twist angle, a: link length, 
            d: link offset, q: joint angle'''
            T = Matrix([[             cos(q),             -sin(q),           0,               a],
                        [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
                        [sin(q) * sin(alpha), cos(q) * sin(alpha),  cos(alpha),  cos(alpha) * d],
                        [                  0,                   0,           0,               1]])
            return T 

        
        def rot_x(q):
            ''' Rotation matrix along x axis'''
            R_x = Matrix([[      1,      0,      0],
                          [      0, cos(q), -sin(q)],
                          [      0, sin(q),  cos(q)]])
            
            return R_x
    
        def rot_y(q):
            ''' Rotation matrix along y axis'''
            R_y = Matrix([[ cos(q),     0, sin(q)],
                          [      0,     1,      0],
                          [-sin(q),     0, cos(q)]])
            
            return R_y

        def rot_z(q):
            ''' Rotation matrix along z axis'''
            R_z = Matrix([[ cos(q),-sin(q), 0],
                          [ sin(q), cos(q), 0],
                          [      0,      0, 1]])
            
            return R_z

        def Euler_angles_from_matrix_URDF_old(R, angles_pre):
            ''' Calculate q4-6 from R3_6 rotation matrix
            Input R is 3x3 rotation matrix, output are Euler angles :q4, q5, q6'''
            r12, r13 = R[0,1], R[0,2]
            r21, r22, r23 = R[1,0], R[1,1], R[1,2] 
            r32, r33 = R[2,1], R[2,2]
            # # Euler angles from rotation matrix
            # q5 = atan2(sqrt(r13**2 + r33**2), r23)
            # q4 = atan2(r33, -r13)
            # q6 = atan2(-r22, r21)
            if np.abs(r23) is not 1:
                q5 = atan2(sqrt(r13**2 + r33**2), r23)
                if sin(q5) < 0:
                    q4 = atan2(-r33, r13)
                    q6 = atan2(r22, -r21)
                else:
                    q4 = atan2(r33, -r13)
                    q6 = atan2(-r22, r21)
            else:
                q6 = angles_pre[5]
                if r23 == 1:
                    q5 = 0
                    q4 = -q6 + atan2(-r12, -r32)
                else:
                    q5 = 0
                    q4 = q6 - atan2(r12, -r32)

            return np.float64(q4), np.float64(q5), np.float64(q6)
        
        def Euler_angles_from_matrix_URDF(R, angles_pre):
            ''' Calculate q4-6 from R3_6 rotation matrix
            Input R is 3x3 rotation matrix, output are Euler angles :q4, q5, q6'''
            r12, r13 = R[0,1], R[0,2]
            r21, r22, r23 = R[1,0], R[1,1], R[1,2] 
            r32, r33 = R[2,1], R[2,2]
            if np.abs(r23) is not 1:
                q5 = np.arctan2(np.sqrt(r13**2 + r33**2), r23)
                if sin(q5) < 0:
                    q4 = np.arctan2(-r33, r13)
                    q6 = np.arctan2(r22, -r21)
                else:
                    q4 = np.arctan2(r33, -r13)
                    q6 = np.arctan2(-r22, r21)
            else:
                q6 = angles_pre[5]
                if r23 == 1:
                    q5 = 0
                    q4 = -q6 + np.arctan2(-r12, -r32)
                else:
                    q5 = 0
                    q4 = q6 - np.arctan2(r12, -r32)

            return q4, q5, q6
        


        def angle_simplify(theta):
            ''' Transfer the theta into [-2pi, 2pi]'''
            return np.float64(np.abs(theta)%(2*np.pi) * np.sign(theta))
        

        # Define variables
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
        r, p, y = symbols("r p y") # end-effector orientation
        px_s, py_s, pz_s = symbols("px_s py_s pz_s") # end-effector position
        R_corr = rot_z(pi) * rot_y(-pi/2) # Compensation matix from URDF to world frame

        theta1,theta2,theta3,theta4,theta5,theta6 = 0,0,0,0,0,0
        angles_pre = (0,0,0,0,0,0)
        r2d = 180./np.pi
        loop_times = []

        # The Modified DH params
        s = {alpha0:       0, a0:      0, d1:  0.75,
             alpha1:   -pi/2, a1:   0.35, d2:     0, q2: q2-pi/2,
             alpha2:       0, a2:   1.25, d3:     0,
             alpha3:   -pi/2, a3: -0.054, d4:   1.5,
             alpha4:    pi/2, a4:      0, d5:     0,
             alpha5:   -pi/2, a5:      0, d6:     0,
             alpha6:       0, a6:      0, d7: 0.303, q7: 0}

        # Define Modified DH Transformation matrix
        if not os.path.exists("R0_3_inv.p"):
            print("-----------------------------------------")
            print("No DH matrices, create and save it.")
            # base_link to link1
            T0_1 = DH_transform(q1, d1, alpha0, a0).subs(s)
            # linke1 to link 2
            T1_2 = DH_transform(q2, d2, alpha1, a1).subs(s)
            # link2 to link3
            T2_3 = DH_transform(q3, d3, alpha2, a2).subs(s)
            # link3 to link4
            T3_4 = DH_transform(q4, d4, alpha3, a3).subs(s)
            # link4 to link5
            T4_5 = DH_transform(q5, d5, alpha4, a4).subs(s)
            # link5 to link6
            T5_6 = DH_transform(q6, d6, alpha5, a5).subs(s)
            # link6 to end-effector
            T6_7 = DH_transform(q7, d7, alpha6, a6).subs(s)
   
            # Create individual transformation matrices
            T0_2 = simplify(T0_1 * T1_2)
            p2_0_sym = T0_2 * Matrix([0,0,0,1])
            # R0_3 inv matrix would be used in R3_6
            T0_3 = simplify(T0_2 * T2_3)
            R0_3 = T0_3[0:3, 0:3]
            R0_3_inv = simplify(R0_3 ** -1)
            # T0_4 = simplify(T0_3 * T3_4)
            # T0_5 = simplify(T0_4 * T4_5)
            # T0_6 = simplify(T0_5 * T5_6)
            # T0_G = simplify(T0_6 * T6_7)
            # R_corr = rot_z(pi) * rot_y(-pi/2)
            #T_total = simplify(T0_G * R_corr)
            # T3_6 = simplify(T3_4*T4_5*T5_6)
            R0_g_sym = simplify(rot_z(y) * rot_y(p) * rot_x(r))
            # pickle.dump(T0_2, open("T0_2.p", "wb"))
            # pickle.dump(T3_6, open("T3_6.p", "wb"))
            pickle.dump(p2_0_sym, open("p2_0_sym.p", "wb"))
            pickle.dump(R0_3_inv, open("R0_3_inv.p", "wb"))
            pickle.dump(R0_g_sym, open("R0_g_sym.p", "wb"))
            
            print("Transformation matrices have been saved!!")
            print("-----------------------------------------")
        else:
            # T0_2 = pickle.load(open("T0_2.p", "rb"))
            # T3_6 = pickle.load(open("T3_6.p", "rb"))
            p2_0_sym= pickle.load(open("p2_0_sym.p", "rb"))
            R0_3_inv = pickle.load(open("R0_3_inv.p", "rb"))
            R0_g_sym = pickle.load(open("R0_g_sym.p", "rb"))
            
            print("-----------------------------------------")
            print("Transformation matrices have been loaded!")         

        # Joint angles calculation
        pg_0 = Matrix([[px_s],[py_s],[pz_s]])
        R0_g = R0_g_sym[0:3,0:3] * R_corr 
        pwc_0 = pg_0 - 0.303 * R0_g * Matrix([[0],[0],[1]])
        theta1_sym = atan2(pwc_0[1], pwc_0[0]).subs(s)

        
        pwc_2 = pwc_0 - p2_0_sym[0:3,:]
        l23 = a2
        l35 = sqrt(a3**2 + d4**2)
        # l25 = sqrt(np.sum(np.square(pwc_2)))
        l25 = sqrt(pwc_2[0]**2 + pwc_2[1]**2 + pwc_2[2]**2)

        # Calculate theta2
        theta22 = atan2(pwc_2[2], sqrt(pwc_2[0]**2 + pwc_2[1]**2))
        c523 = (-l35**2 + l23**2 + l25**2) / (2*l23 * l25)
        theta21 = atan2(sqrt(1 - c523**2), c523)
        theta2_sym = (pi/2 - (theta21 + theta22)).subs(s)

        # Calculate theta3
        theta31 = atan2(a3,d4)
        c235 = (l25**2 - l23**2 - l35**2) / (2*l23*l35)
        theta32 = atan2(sqrt(1-c235**2), c235)
        theta3_sym = (theta32 + theta31 - pi/2).subs(s)

        # Transfer symbol into numpy function to evaluate an expression more efficient
        theta1_f = lambdify((px_s,py_s,pz_s,r,p,y), theta1_sym, "numpy")
        theta2_f = lambdify((q1,px_s,py_s,pz_s,r,p,y), theta2_sym, "numpy")
        theta3_f = lambdify((q1,px_s,py_s,pz_s,r,p,y), theta3_sym, "numpy")

        # Calculate R3_6 for theta4-6 
        R3_6_sym = R0_3_inv * R0_g
        R3_6_f = lambdify((q1,q2,q3,px_s,py_s,pz_s,r,p,y), R3_6_sym, "numpy")


        loop_start_time = time.time()
        for x in xrange(0, len(req.poses)):
            loop_current_time = time.time()
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()
            
            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            # Calculate joint angles using Geometric IK method           
            theta1 = theta1_f(px,py,pz,roll, pitch, yaw)          
            theta2 = theta2_f(theta1, px,py,pz,roll, pitch, yaw)
            theta3 = theta3_f(theta1, px,py,pz,roll, pitch, yaw)

            R3_6 = R3_6_f(theta1, theta2, theta3,px,py,pz,roll, pitch, yaw)
            theta4, theta5, theta6 = Euler_angles_from_matrix_URDF(R3_6, angles_pre)

            angles_pre = (theta1, theta2, theta3, theta4, theta5, theta6)

            loop_time = time.time()-loop_current_time
            loop_times.append(loop_time)
            # print("Calculating trajectory:{:2d}, Time cost:{:>2f}".format(x, loop_time))

            # Populate response for the IK request
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
        print("Inverse kinematics calculation has been done!")
        print("Total calculation time: {:>4f}".format(time.time() - loop_start_time))
        print("Average calculation time: {:>4f}".format(np.mean(loop_times)))
        print("-----------------------------------------")
        print("theta1: {:>4f}".format(theta1 * r2d))
        print("theta2: {:>4f}".format(theta2 * r2d))
        print("theta3: {:>4f}".format(theta3 * r2d))
        print("theta4: {:>4f}".format(theta4 * r2d))
        print("theta5: {:>4f}".format(theta5 * r2d))
        print("theta6: {:>4f}".format(theta6 * r2d))
        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        print("-----------------------------------------")
        print("Moving robot arm...")
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print ("Ready to receive an IK request")
    rospy.spin()

if __name__ == "__main__":
    IK_server()
