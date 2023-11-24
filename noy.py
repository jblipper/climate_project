from math import sin,cos,atan2,pi
point_sel=(-0.118092,51.509865)
def get_D(point):
     R=6371000 # radius of Earth(m)
     x1=point_sel[0]*pi/180
     y1=point_sel[1]*pi/180
     x2=point[0]*pi/180
     y2=point[1]*pi/180
     a=sin((y2-y1)/2)**2+cos(y1)*cos(y2)*sin((x2-x1)/2)**2
     c=2*atan2(a**0.5,(1-a)**0.5)
     d=R*c
     return d
point=(-114.066666,51.049999)
print(get_D(point)/1000)
