# %%
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 13:54:32 2024
#1 can probably take more than an hour to finish 3 for loops in the force function!!
@author: Hins and his lovely companions
"""# importing library and intial values
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
#from numpy import *
import random
import os

L=64
#L=32
N=L**2
restlen=np.ones(shape=((L**2),6))
p=0.8
bond_len = 1
dipole_number=2
#r=6
r=25
#center=120
#center=496
center=2080
final_restlen=0.9
total_step=10000
#total_step=1    
base_directory = f'E:/FIRE_algorithm/Test_run/max steps={total_step}_circle{L}/{dipole_number}dipoles/p={p}'
#base_directory = '/home/abhinav/david/Hins/test_run/'
if not os.path.exists(base_directory):
    os.makedirs(base_directory)
# %%
#******************************************************************************************************
#the connection matrix for each node
def connect_matrix(L):
    connect=-1*np.ones((L**2,6), dtype=int)
    n=0
    for n in range(L**2):
        connect[n][0]=n-1
        connect[n][1]=n-L-1
        connect[n][2]=n-L
        connect[n][3]=n+1
        connect[n][4]=n+L
        connect[n][5]=n+L-1
    #left edges even rows
        if n%L==0:
            connect[n][0]=-1
    #even rows        
        if (int(n/L)%2)==0:
            connect[n][1]=n-L
            connect[n][2]=n-L+1
            connect[n][4]=n+L+1
            connect[n][5]=n+L
            #bottom nodes      
            if n<L:
                connect[n][1]=-1
                connect[n][2]=-1
            #top nodes
        if (n+L) >= (L**2):
            connect[n][4]=-1
            connect[n][5]=-1
    #left edges, odd rows        
        if (((n/L)-1))%2==0:
            connect[n][1]=-1
            connect[n][5]=-1
    #right edges odd rows
        if (n+1)%L==0:
            connect[n][3]=-1
    #right edges, even rows
        if (((n+1)/L)-1)%2==0:
            connect[n][2]=-1
            connect[n][4]=-1
    #connect_full=connect.copy()
    return  connect
connect=connect_matrix(L)  
#************************************************************************************************************
 # initialize the network at some positions xcord and ycord
def initial_network(L):
    xcord=np.zeros(shape=(L**2))
    for i in range(N):
        #for r in range(L):
            if (i//L)%2==0:
                xcord[i] =i%L+bond_len*1.5
            else:
                xcord[i]= i%L+bond_len
    #j=[]
    initial_value = 1
    ycord=np.zeros(shape=(L**2))
    for j in range(N):
        if j<L:
            ycord[j]=initial_value
        else:
            ycord[j] = initial_value +(j//L)*0.5*np.sqrt(3)
    return xcord, ycord

xcord,ycord=initial_network(L)

#***********************************************************************************
def circle_boundary(xcord,ycord):
    circle_boundary_node=[]
    inner_circle=[]
    connect_circle=connect.copy()
    radius=[]
    #outer_circle=[]
    x_radius=xcord-xcord[center]
    y_radius=ycord-ycord[center]
    radius=(np.sqrt(x_radius**2 + y_radius**2))
    for i in range(N):
            if radius[i]<=r-0.5:
              inner_circle.append(i)  
            #if (r-1)<=radius<=r:
            if (r-0.5)<=radius[i]<=r+0.5:
              circle_boundary_node.append(i)
            for j in range(6):
                    if radius[connect_circle[i][j]]>r+0.5 or radius[i]>r+0.5:
                        #outer_circle.append(i)
                        connect_circle[i][j]=-1
    return radius,inner_circle, circle_boundary_node,connect_circle
radius,inner_circle, circle_boundary_node,connect_circle=circle_boundary(xcord,ycord)
# %%
#*****************************************************************************************************

dip_center=[]
dipole_circle=[]
def create_dipoles(N): 
    for i in range(N):
        if radius[i]<=r-12.5:
          dipole_circle.append(i)
    #rand_dipole=(random.sample(inner_circle,3))
    while len(dip_center) < dipole_number:
        new_dipole = np.random.choice(dipole_circle, size=1)[0]
        
        # Check if the new number is far enough from all existing numbers
        if all(abs(new_dipole - num) >= L for num in dip_center):
            dip_center.append(new_dipole)
    #dip_center =(random.sample(dipole_circle,dipole_number))
   #dip_center =(np.random.choice(inner_circle,dipole_number,replace=False))
    dip_node=dip_center.copy()
    for m in range(6):
      for c in range(0,len(dip_center)):
             dip_node.append(connect[dip_center[c]][m])
    inner_bond = [num for num in inner_circle if num not in dip_center]
    return dip_center, dip_node,inner_bond
dip_center, dip_node, inner_bond=create_dipoles(N) 
# %%
#*****************************************************************************************************

# Trying to keep the boundaries and dipoles out of the removal!!!
## Function to generate a random combination
def generate_random_combination(min_val2, max_val2):
    rand_node = np.random.choice(inner_bond, size=1)[0]
    rand_bond = random.randint(min_val2, max_val2)
    return rand_node, rand_bond

# val1 corresponds to random node number being picked;
# val2 corresponds to random bond number being removed
#min_val1=285
#max_val1=L**2-L
min_val2=0
max_val2=5


# Main function to generate unique bond removal
def random_bond_removal(connect_circle, p):
    bond_count=int((L**2*6-np.sum(connect_circle==-1))/2)
    remove_bonds=int(bond_count*(1-p))
    unique_combinations = []
    while len(unique_combinations) < remove_bonds:
        node_index, bond_index = generate_random_combination(min_val2, max_val2)
        combination=(node_index, bond_index)
        if connect_circle[node_index][bond_index] not in circle_boundary_node and connect_circle[node_index][bond_index] not in dip_center:
            #print(node_index)
            unique_combinations.append(combination)
            othernode_index=connect_circle[node_index][bond_index]
            #print(combination)
            connect_circle[node_index][bond_index]=-1
            for z in range(0,6):
                if connect_circle[othernode_index][z]==node_index:
                    connect_circle[othernode_index][z]=-1
    return connect_circle
connect_circle=random_bond_removal(connect_circle, p)

#**********************************************************************************************
# plotting function below
def network_plot(xnew,ynew,key):
        #threshold = 1e-4
        #threshold = 1e-5
        threshold = 1e-6
        #threshold = 1e-7
        #threshold = 1e-8
        #threshold = 1e-9
        linesvec=[]
        BLenvec=[]
        colorsvec=[]
        #circle=plt.Circle(plt.Circle((0.5, 0.5), 0.3, color='blue', fill=False))
        # change the following lines variable and make it suitable for your lattice of nodes.
        for n in range(0,L**2):
                for m in range (0,6):
                     if connect_circle[n][m]!=-1:
                        x1=xnew[connect_circle[n][m]]
                        x2=xnew[n]
                        y1=ynew[connect_circle[n][m]]    
                        y2=ynew[n]
                        lines = [(x1, y1), (x2, y2)]
                        BLen=(((x2-x1)**2+(y2-y1)**2)**0.5)
                        BLenvec.append(BLen)
                        linesvec.append(lines)
        # format of lines: [[x1,y1],[x2,y2]]
                        if BLen-restlen[n][m]>threshold:
                            colors = 'b'
            #            elif BLenvec[l]-1<0:
                        elif BLen-restlen[n][m]<-1*threshold:
                            colors = 'r'
                        else:
                            colors = 'k'
                        colorsvec.append(colors)
                    #colors = 'k'
        lc = LineCollection(linesvec,linewidths=0.5, colors=colorsvec)
        fig = plt.figure(dpi = 300)
        #fig = plt.figure(dpi=1000)
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.add_collection(lc)
        #ax1.add_patch(circle)
        oneD_x=xnew.flatten()
        oneD_y=ynew.flatten()
        color_sc = np.full((L**2), 'grey')
        color_sc[circle_boundary_node] = 'c'
        color_sc[dip_node] = 'g'
        #sc_dipole_number = np.full((L**2),20)
        sc_dipole_number = np.full((L**2),1)
        #sc_dipole_number[dip_node] =70
        sc_dipole_number[dip_node] =15
        sc_dipole_number[circle_boundary_node] =10
        #print(oneD_x)
        # change the scatter plot. Give it x and y locations of each node in array format.
        ax1.scatter([oneD_x],[oneD_y], color = color_sc, s=sc_dipole_number)

        ax1.autoscale()
        ax1.set_title('Network')
        # plt.show()
        
        plt.savefig(base_directory+f"{dipole_number} network_"+key+f"step{total_step}"+".png")

network_plot(xcord,ycord, 'initial')
# %%


#***********************************************************************************
# creating dx and dy matrix: for every connected node they include the distance from a node
# using the distances to find forces from every neighbor and then sum all the x and y components of the force.
def  restlen_change(restlen):
    # final_restlen=0.8
    for n in range (0,L**2):
        for m in range (0,6):
            for c in range(0,len(dip_center)):
                if connect_circle[n][m]==dip_center[c]:
                    restlen[n][m]=final_restlen
    #Drestlen[dip_center,:]=0.01 
    restlen[dip_center,:]=final_restlen      
    #if restlen[dip_center[0]][0]>final_restlen:
     #   restlen-=Drestlen
    return restlen
restlen=restlen_change(restlen)
#the new length of each spring
xnew = xcord.copy()
ynew = ycord.copy()
def force(xnew,ynew):
    #restlen=restlen_change(restlen)
    xspring_len=np.zeros(shape=(L**2,6))
    yspring_len=np.zeros(shape=(L**2,6))
    dspring=np.zeros(shape=(L**2,6))
    Fx_spring=np.zeros(shape=((L**2),6)); 
    Fy_spring=np.zeros(shape=((L**2),6));
    #radialF_boundary=np.zeros(shape=len(inner_circle))
    #dip_node=[120,121,135,136,137,152,153]
    #dip_node=[135, 136, 137]
    #for i in range (N):
    for i in range(N):
      for j in range (6):
        #if connect[inner_circle[c]][j]!=-1:
        if connect_circle[i][j]!=-1:
           # below is dx = x_m - x_n i.e. the difference of node connect_circleed by mth spring - nth node
           xspring_len[i][j] = xnew[connect_circle[i][j]] - xnew[i]
           yspring_len[i][j] = ynew[connect_circle[i][j]] - ynew[i]
           spring_len=np.sqrt(np.square(xspring_len)+np.square(yspring_len))
        #if radius[connect_circle[[i]][j]]<r+0.5:
           dspring[i][j]=spring_len[i][j]-restlen[i][j]
           Fx_spring[i][j]=(k*dspring[i][j]*(xspring_len[i][j]/spring_len[i][j]))
           Fy_spring[i][j]=(k*dspring[i][j]*(yspring_len[i][j]/spring_len[i][j])) 
        #if radius[connect_circle[inner_circle[c]][j]]>r+0.5:
         #   Fx_spring[inner_circle[c]][j]=0
          #  Fy_spring[inner_circle[c]][j]=0
        
    Fx=Fx_spring.sum(axis=1)
    Fy=Fy_spring.sum(axis=1)
    #Fx = np.reshape(sumFx,(L,L))
    #Fy = np.reshape(sumFy,(L,L))
    PE=np.sum(1/4*k*np.square(dspring))
    #PE=sum(PE)
    #PE=PE.sum(axis=0)
    sum_F=np.sqrt(np.sum(np.square(Fx)+np.square(Fy))) 
    return Fx,Fy,PE,sum_F,spring_len
#**********************************************************************************
def molecular_dynamics(xnew,ynew,v_x,v_y):
    a_x=np.zeros(shape=(N))
    a_y=np.zeros(shape=(N))
    for n in range (N):
        #excluding boundary nodes
     # if any(circle_boundary_node)==n:
        
        if radius[n]<(r-0.5):
            a_x[n]=Fx[n]/m; a_y[n]=Fy[n]/m;
            xnew[n]+=v_x[n]*dt+1/2*a_x[n]*(dt**2);
            ynew[n]+=v_y[n]*dt+1/2*a_y[n]*(dt**2);
            v_x[n]+=1/2*dt*(a_x[n]+(Fx[n]/m));
            v_y[n]+=1/2*dt*(a_y[n]+(Fy[n]/m));
        else:
            xnew[n] = xcord[n]
            ynew[n] = ycord[n]
            v_x[n] = 0
            v_y[n] = 0

    return xnew,ynew,v_x,v_y
#**********************************************************************************
#FIRE constants
N_min=5
f_inc=1.1
f_dec=0.5
a_start=0.1
f_a=0.99
dt=0.001
dt_max=10*dt
m=1
k=0.1

# initial matrix setup
v_x=np.zeros(shape=(N))
v_y=np.zeros(shape=(N))
a_xvec=np.zeros(shape=((L**2),6)); 
a_yvec=np.zeros(shape=((L**2),6));
v_xvec=[]; 
v_yvec=[];
sum_velocity=np.transpose([v_x, v_y]).reshape(L**2,2)
P_vec=np.zeros(shape=((L**2),6)); PE_vec=[]
#different horizon(network)
cut=0
steps=0
a=a_start

#FIRE
axvec = []
ayvec = []
xvec = []
yvec = []
F_vec=[]
sumF_vec=[]
AvgF_vec=[]
list_total_p = []
F_threshold_per_node=2.5e-7

xnew_inner = xcord[inner_circle]
ynew_inner = ycord[inner_circle]

while steps<total_step:
    #if radius<=r:
        #restlen=restlen_change(restlen)
        if steps%1000 == 0:
            print('step number = ', steps)
        Fx,Fy,PE,sum_F,spring_len=force(xnew,ynew)
        Avg_sumF=sum_F/(np.shape(inner_circle))
        sumF_vec=np.append(sumF_vec,sum_F)
        AvgF_vec=np.append(AvgF_vec,Avg_sumF)
        PE_vec=np.append(PE_vec,PE)
        V=np.sqrt((np.square(v_x)+np.square(v_y)))
        P=v_x*Fx+v_y*Fy
        P_vec = np.append(P_vec,P)
        Total_P=np.sum(P)
        list_total_p.append(Total_P)
        #if Avg_sumF<F_threshold_per_node:
          #  break
        if Total_P<=0:
            print("P is ZERO!! Step number = ", steps)
            cut=steps
            #print(cut)
            v_x=np.zeros(shape=(N))
            v_y=np.zeros(shape=(N))
            dt*=f_dec
            a=a_start
        elif Total_P>0 and (steps-cut)>N_min:
            #for n in range (0,L**2):
              #if n<L or (n+L)>=(L**2) or n%L==0 or (((n/L)-1))%2==0 or (n+1)%L==0 or (((n+1)/L)-1)%2==0:
                  #v_x[n]=0
                  #v_y[n]=0
              #else:
            vx_term1 = (1-a)*v_x
            vx_term2_numerator = a*Fx*V
            vx_term2_denominator = (np.sqrt((np.square(Fx)+np.square(Fy))))
            vx_term2 = np.divide(vx_term2_numerator, vx_term2_denominator, out=np.zeros_like(vx_term2_numerator), where=vx_term2_denominator!=0)
            v_x = vx_term1 + vx_term2
        
            vy_term1 = (1-a)*v_y
            vy_term2_numerator = a*Fy*V
            vy_term2_denominator = (np.sqrt((np.square(Fx)+np.square(Fy))))
            vy_term2 = np.divide(vy_term2_numerator, vy_term2_denominator, out=np.zeros_like(vy_term2_numerator), where=vy_term2_denominator!=0)
            v_y = vy_term1 + vy_term2
            dt = min( dt*f_inc, dt_max )
            a = a * f_dec
    #Molecular Dynamics(MD)
        xnew,ynew,v_x,v_y=molecular_dynamics(xnew, ynew, v_x, v_y)
        v_xvec = np.append(v_xvec,v_x); v_yvec = np.append(v_yvec,v_y)
        #axvec = np.append(axvec,a_x); ayvec = np.append(ayvec,a_y)
        xvec = np.append(xvec,xnew); yvec = np.append(yvec,ynew)
        steps+=1
#xvec_new=np.reshape(xvec,(i,L**2))
#yvec_new=np.reshape(yvec,(i,L**2))
#Fx_Fy=np.column_stack((Fx,Fy))
#xnew_ynew=np.column_stack((xnew,ynew))
file_Fx = "Fx.npy"
file_Fy = "Fy.npy"
#file_Fx_Fy="Fx_Fy.npy"
#file_xnew_ynew="xnew_ynew.npy"
file_sumF="sumF.npy"
file_xnew = "xnew_.npy"
file_ynew = "ynew_.npy"
file_xinit = "x_init.npy"
file_yinit = "y_init.npy"
file_restlen= "restlen.npy"
file_PE_vec = "PE_vec.npy"
file_AvgF_vec = "AvgF_vec.npy"
file_innercircle= "innercircle.npy"
file_boundarynode = "circle_boundary_node.npy"
file_springlen="springlen.npy"
file_dipnode="dipole_node.npy"
file_dipcenter="dipole_center.npy"

file_path_Fx = os.path.join(base_directory, file_Fx)
file_path_Fy = os.path.join(base_directory, file_Fy)
file_path_xnew = os.path.join(base_directory, file_xnew)
file_path_ynew = os.path.join(base_directory, file_ynew)
#file_path_Fx_Fy = os.path.join(base_directory, file_Fx_Fy)
#file_path_xnew_ynew = os.path.join(base_directory, file_xnew_ynew)
file_path_restlen = os.path.join(base_directory, file_restlen)
file_path_sumF = os.path.join(base_directory, file_sumF)
file_path_PE = os.path.join(base_directory, file_PE_vec)
file_path_xinit = os.path.join(base_directory, file_xinit)
file_path_yinit = os.path.join(base_directory, file_yinit)
file_path_AvgF = os.path.join(base_directory, file_AvgF_vec)
file_path_innercircle = os.path.join(base_directory, file_innercircle)
file_path_boundarynode = os.path.join(base_directory, file_boundarynode)
file_path_springlen = os.path.join(base_directory, file_springlen)
file_path_dipnode = os.path.join(base_directory, file_dipnode)
file_path_dipcenter = os.path.join(base_directory, file_dipcenter)

np.savetxt(file_path_Fx,Fx, delimiter=',')
np.savetxt(file_path_Fy,Fy, delimiter=',')
np.savetxt(file_path_xnew, xnew, delimiter=',') 
np.savetxt(file_path_ynew, ynew, delimiter=',')
#np.savetxt(file_path_Fx_Fy,Fx_Fy, delimiter=',')
#np.savetxt(file_path_xnew_ynew,xnew_ynew, delimiter=',')
np.savetxt(file_path_xinit, xcord, delimiter=',') 
np.savetxt(file_path_yinit, ycord, delimiter=',')
np.savetxt(file_path_restlen, restlen, delimiter=',')
np.savetxt(file_path_PE, PE_vec, delimiter=',')
np.savetxt(file_path_AvgF, AvgF_vec, delimiter=',')
np.savetxt(file_path_sumF, np.sqrt(Fx**2+Fy**2), delimiter=',')
np.savetxt(file_path_innercircle, inner_circle, delimiter=',')
np.savetxt(file_path_boundarynode, circle_boundary_node, delimiter=',')
np.savetxt(file_path_springlen, spring_len, delimiter=',')
np.savetxt(file_path_dipnode, dip_node, delimiter=',')
np.savetxt(file_path_dipcenter, dip_center, delimiter=',')

print("Wohoooo it ran!!")
#**********************************************************************************
#call for plots
plt.figure()
plt.plot(PE_vec)
plt.title('potential energy')
plt.savefig(base_directory+ f"potential energy step{total_step}.png")
#plt.show()
print(f'final potential energy value={PE_vec[-1]}')
print(PE_vec[-2])
print(((PE_vec[-1]-PE_vec[-2])/PE_vec[-1])*100)
#print("shape of potential energy = ", np.shape(PE_vec))
plt.figure()
plt.plot(sumF_vec)
plt.title('force of the network')
plt.savefig(base_directory+ f"network force step{total_step}.png")
#plt.show()
#plt.plot(PE_vec[305:330])
#print(PE_vec[240:260])
xaxis=np.arange(0,2000)
plt.figure()
plt.plot(list_total_p)
plt.title('power of the network')
plt.savefig(base_directory+ f"power step {total_step}.png")
#plt.show()
print(min(list_total_p))
print(np.shape(np.where(np.array(list_total_p)<=0)))
print(np.shape(list_total_p))

network_plot(xnew,ynew, 'final')

comparison=xcord[circle_boundary_node]==xnew[circle_boundary_node]
print(comparison)
F_net=np.sqrt(np.square(Fx)+np.square(Fy))
print(F_net)
# %%

'''
#**********************************************************************************
# calculating the diople moment at the boundary
# input_directory = '/home/abhinav/david/Hins/1000steps_circular/'
# input_directory = '/home/abhinav/david/Hins/test_run/'
#input_directory = '/home/abhinav/david/Hins/test_run/8000steps/'
input_directory = f'E:/FIRE_algorithm/Test_run2/max steps={total_step}_test_circle{L}/2dipoles/'
fx_file = input_directory + 'Fx.npy'
fy_file = input_directory + 'Fy.npy'

xpos_file = input_directory + 'xnew_.npy'
ypos_file = input_directory + 'ynew_.npy'

fx_data = np.loadtxt(fx_file)
fy_data = np.loadtxt(fy_file)
xpos_data = np.loadtxt(xpos_file)
ypos_data = np.loadtxt(ypos_file)

fx_boundary = fx_data[circle_boundary_node]
fy_boundary = fy_data[circle_boundary_node]
xpos_boundary = xpos_data[circle_boundary_node]
ypos_boundary = ypos_data[circle_boundary_node]

center_x = xpos_data[center]
center_y = ypos_data[center]

radial_vec = []
for i in range(0,np.shape(xpos_boundary)[0]):
    radial_x = center_x - xpos_boundary[i]
    radial_y = center_y - ypos_boundary[i]
        
    radial_vec.append([radial_x,radial_y])

radial_vec = np.array(radial_vec)

bndry_force_vec = np.column_stack((fx_boundary, fy_boundary))

xx_moment_sum_array = bndry_force_vec[:,0] * radial_vec[:,0]    # fx dot rx
yy_moment_sum_array = bndry_force_vec[:,1] * radial_vec[:,1]    # fy dot ry

xx_moment_sum = np.sum(xx_moment_sum_array)
yy_moment_sum = np.sum(yy_moment_sum_array)

dfar_fire = xx_moment_sum + yy_moment_sum     # this is the dipole moment at the boundary

#**********************************************************************************

# d_loc calculations:: it is the local dipole moment

import math
mu = 1
mu_c = 1
def angle(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2-y1
    return(math.atan2(dy,dx))

# validating that the above short-cut is correct by finding unbalanced forces on the outer dipole nodes
y_force_unbalanced = []
x_force_unbalanced = []
strain_loc = []
test_force_list = []
dist_dip_node_list = []
bonds_list = [0,1,2,3,4,5]

j = 0
for i in range(0,len(dip_node)):
    if i%7 == 0: 
        last_center = dip_node[i]
        index = int(i/7)
        y_force_unbalanced.append([])
        x_force_unbalanced.append([])
        dip_nodes_outer = dip_node[i+1:i+7]
        test_force_list.append([])
        
        
        for bond in bonds_list:
            conn_val = connect[dip_node[i]][bond]     # node that is the fourth bond connnected to
            x1 = xpos_data[dip_node[i]]
            y1 = ypos_data[dip_node[i]]
            x2 = xpos_data[conn_val]
            y2 = ypos_data[conn_val]
            dx = x1-x2
            dy = y1-y2
            dist = np.sqrt(dx**2+dy**2) 
            dist_dip_node_list.append(dist)
            print('dist is: ', dist)

                
    if i%7 != 0:

        # find the orientation of the spring connected to dipole node
        # x_last_center = xpos[-1,last_center]
        # y_last_center = ypos[-1,last_center]
        x_last_center = center_x
        y_last_center = center_y
        # x_node = xpos[-1,dip_nodes[i]]
        # y_node = ypos[-1,dip_nodes[i]]
        x_node = xpos_data[dip_node[i]]
        y_node = ypos_data[dip_node[i]]
        
        spring_vec = [x_node-x_last_center, y_node-y_last_center]
        
        angle_spring = angle(x_node, y_node, x_last_center, y_last_center)
        
        force_node_y = 0
        force_node_x = 0
        test_force_node = 0
        
        for bond in bonds_list:
            # strain_val = strain[-1][dip_node[i]][bond]    # fourth bond strain
            # conn_val = conn_node[dip_node[i]][bond]     # node that is the fourth bond connnected to
            conn_val = connect[dip_node[i]][bond]     # node that is the fourth bond connnected to
            
            # if np.isnan(strain_val) == False:# and conn_val != last_center:
            if conn_val != -1: # and conn_val != last_center:
                
                # print(dip_nodes[i], conn_val, strain_active_force)
                
                x1 = xpos_data[dip_node[i]]
                y1 = ypos_data[dip_node[i]]
                x2 = xpos_data[conn_val]
                y2 = ypos_data[conn_val]
                dx = x1-x2
                dy = y1-y2
                dist = np.sqrt(dx**2+dy**2) 
                strain_val = dist - restlen[dip_node[i],bond]
                angle_val = angle(x1, y1, x2, y2)

            if conn_val == last_center:
                # strain_active_force = 1. - (rlen + strain_val)
                strain_active_force = 1. - (final_restlen + strain_val)
                # print(dip_nodes[i], bond, strain_active_force, strain_val)
            else:
                strain_active_force = strain_val

                if strain_val >= 0:
#                    force_bond = -1*strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
#                    force_bond = strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                    force_bond = strain_active_force*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                elif strain_val < 0:
#                    force_bond = -1*strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
#                    force_bond = strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
                    force_bond = strain_active_force*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary

                force_bond_y = force_bond*np.sin(angle_val)
                force_bond_x = force_bond*np.cos(angle_val)
                
                # if dip_nodes[i] == 2079:
                #     print(bond, force_bond_x, force_bond_y, np.sqrt(np.square(force_bond_x) + np.square(force_bond_y)), force_bond)
                
                force_node_y = force_node_y + force_bond_y  # adding bond forces to force on nodes
                force_node_x = force_node_x + force_bond_x  # adding bond forces to force on nodes
                
                force_unbalanced_temp = force_bond*np.cos(angle_spring)
                # finding the projection of spring force along the radial outward vector                
                test_force = np.dot([force_bond_x, force_bond_y], spring_vec) / np.linalg.norm(spring_vec)
                # test_force = np.dot([force_bond_x, force_bond_y], spring_vec)
                # test_force = 
                
                if i%7 == 1 and bond == 3:
                    test_force = -1*test_force
                if i%7 == 2 and bond == 4:
                    test_force = -1*test_force
                if i%7 == 3 and bond == 5:
                    test_force = -1*test_force
                if i%7 == 4 and bond == 0:
                    test_force = -1*test_force
                if i%7 == 5 and bond == 1:
                    test_force = -1*test_force
                if i%7 == 6 and bond == 2:
                    test_force = -1*test_force
                
                if dip_node[i] == 2079:
                    print([force_bond_x, force_bond_y])#, spring_vec)
                    # print(bond, force_bond, angle_val, force_unbalanced_temp, test_force, np.dot([force_bond_x, force_bond_y], spring_vec), np.linalg.norm(spring_vec))
                    # print(test_force)
                
                test_force_node = test_force_node + test_force
                    
        test_force_list[index].append(test_force_node)
                                
        # y_force_unbalanced[index].append(force_node_y)
        # x_force_unbalanced[index].append(force_node_x)
        
#        y_force_unbalanced[i].append(force_node_y)
#        x_force_unbalanced[i].append(force_node_x)

        j += 1

test_force_list = np.array(test_force_list)

# y_force_unbalanced = np.array(y_force_unbalanced)
# x_force_unbalanced = np.array(x_force_unbalanced)

# force_mag_unbalanced_list = np.sqrt(np.square(x_force_unbalanced) + np.square(y_force_unbalanced))
# force_mag_unbalanced = np.sum(force_mag_unbalanced_list)

# finding dipole moment using unblanced forces on the outer nodes
# d_loc_unbalanced_list = dist_dip_node_list*force_mag_unbalanced_list
# d_loc_unbalanced = np.sum(d_loc_unbalanced_list)
force_list=np.reshape(test_force_list,-1)
d_loc_test_array = dist_dip_node_list*force_list
d_loc_test = np.sum(d_loc_test_array)

print('D_loc is : ', d_loc_test)
print('D_far is : ', dfar_fire)
file_dfar= "dfar.npy"
file_d_loc_test = "d_loc_test.npy"

file_path_dfar = os.path.join(base_directory, file_dfar)
file_path_d_loc_test = os.path.join(base_directory, file_d_loc_test)

# %%
np.savetxt(file_path_dfar, [dfar_fire])

np.savetxt(file_path_d_loc_test, [d_loc_test])

dist_dip_node_list
Out[5]: 
[[0.9610280144914007,
  0.9610280144914022,
  0.9610280144914022,
  0.9610280144914043,
  0.9610280144913986,
  0.9610280144913973]]
'''


