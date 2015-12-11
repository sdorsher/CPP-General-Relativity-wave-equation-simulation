import re, os

numelem=5
factor = 2
for count in range(1,10):
    print str(numelem)
    file = open('params.tmpl','r')
    file2 = open('params.cfg','w')
    pnumelem = re.compile('numelems(\s)*=(\s)*(\d)*;')
    numelemstring = "numelems = " + str(numelem) + ";"
    for line in file:
#        print pnumelem.search(line)
        line2 = pnumelem.sub(numelemstring,line)
        file2.write(line2)
    file.close()
    file2.close()
    cmd = '/home/sdorsher/DG1D-CPP-development/dg1D'
    os.system(cmd)
    numelem=numelem*factor
