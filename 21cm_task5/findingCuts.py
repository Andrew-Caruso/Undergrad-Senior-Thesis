import numpy 

L = 1000 #physical distance in Mpc/h
N = 1024 #number of cells corresponding to L

print("Allowed integer cuts:")
#go from 0 to 200 
for i in range(101):
	if i != 0:
		a = L/i
		b = N/i
		a_bool = a.is_integer() 
		b_bool = b.is_integer() 
		if a_bool == 1 and  b_bool == 1:
			print("a_bool, b_bool:",a_bool, " ", b_bool) 
			print("i: ",i)
			print("a:",a," b:",b,"\n")
