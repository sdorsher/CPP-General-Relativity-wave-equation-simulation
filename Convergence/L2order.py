import re, os

for order in range(2,12):
    file = open('params.tmpl','r')
    file2 = open('params.cfg','w')
    porder = re.compile('elemorder(\s)*=(\s)*(\d)*;')
    orderstring = "elemorder = " + str(order) + ";"
    for line in file:
        line2 = porder.sub(orderstring,line)
        file2.write(line2)
    file.close()
    file2.close()
    cmd = '/home/sdorsher/DG1D-CPP-development/dg1D'
    os.system(cmd)
