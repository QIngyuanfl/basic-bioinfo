#! /usr/bin/env python
import sys
import re
print 'Usage <infile> <outfile>'
file = open(sys.argv[1],'r')
out=open(sys.argv[2],'w')
a=[]
dist={}
for i in file:
	if len(i.strip()) != 0:
		if i.strip()[0] == "=":
			if a == []:
				pass
			else:
				a[0]=a[0].split(" ")[1]
				num=re.findall(r'\d+',a[1])
				dist[int(num[0])]=a
			a=[]
		else:
			a.append("%s" % ( i))
	if len(i.strip()) == 0:
		a.append("%s" % ( i))
a[0]=a[0].split(" ")[1]
num=re.findall(r'\d+',a[1])
dist[int(num[0])]=a

c=dist.keys()
y=1
while 1:
	try:
		d=min(c)
		dist[d][0]=str(y)+". "+dist[d][0]
		out.write("%s\n%s" % ("========================================",''.join(dist[d])))
		c.remove(d)
		y=y+1
	except ValueError:
		break
