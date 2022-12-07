import sys
num = 0
if len(sys.argv) >= 2:
   num = int(sys.argv[1])

temp_file = open("temp.txt",'r')
temp = temp_file.readlines()
temp_list = temp[0].split(',')
for i in range(num):
    mdp_pre = open("npt-nopr-md.mdp",'r')
    mdp = open("remd%d/remd.mdp"%(i), 'w')
    lines = mdp_pre.readlines()
    for line in lines:
        newline = line.replace("310", str(temp_list[i]))
        mdp.write(newline)
