import numpy as np

from Law_Of_Motion import LOM_Force

class OctTree():
    def __init__(self,objects,a,b,c,d,e,f):
        self.objects=objects
        self.n_objects=len(self.objects)
        self.children=[None,None,None,None,None,None,None,None]
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.e=e
        self.f=f
        self.m_center()
        self.build()

    def build(self, volume_tolerance=1e-10):
        
        if self.n_objects>1:
            TNW=[]
            TNE=[]
            TSW=[]
            TSE=[]
            BNW=[]
            BNE=[]
            BSW=[]
            BSE=[]
            for i in range(self.n_objects):
                m1=(self.a+self.b)/2
                m2=(self.c+self.d)/2
                m3=(self.e+self.f)/2
                if self.objects[i].position[0]>m1:
                    if self.objects[i].position[1]>m2:
                        if self.objects[i].position[2]>m3:
                            TNE.append(self.objects[i])
                        else:
                            BNE.append(self.objects[i])
                    else :
                        if self.objects[i].position[2]>m3:
                            TSE.append(self.objects[i])
                        else:
                            BSE.append(self.objects[i])
                else :
                    if self.objects[i].position[1]>m2:
                        if self.objects[i].position[2]>m3:
                            TNW.append(self.objects[i])
                        else:
                            BNW.append(self.objects[i])
                    else :
                        if self.objects[i].position[2]>m3:
                            TSW.append(self.objects[i])
                        else:
                            BSW.append(self.objects[i])
                        
            self.children[0]=OctTree(TNW,m1,self.b,self.c,m2,self.e,m3)
            self.children[1]=OctTree(TNE,self.a,m1,self.c,m2,self.e,m3)
            self.children[2]=OctTree(TSW,m1,self.b,m2,self.d,self.e,m3)
            self.children[3]=OctTree(TSE,self.a,m1,m2,self.d,self.e,m3)
            self.children[4]=OctTree(BNW,m1,self.b,self.c,m2,m3,self.f)
            self.children[5]=OctTree(BNE,self.a,m1,self.c,m2,m3,self.f)
            self.children[6]=OctTree(BSW,m1,self.b,m2,self.d,m3,self.f)
            self.children[7]=OctTree(BSE,self.a,m1,m2,self.d,m3,self.f)

    def m_center(self):
        if self.n_objects>0:
            xc=0
            yc=0
            zc=0
            m=0
            for i in range(self.n_objects):
                xc=xc+self.objects[i].mass*self.objects[i].position[0]
                yc=yc+self.objects[i].mass*self.objects[i].position[1]
                zc=zc+self.objects[i].mass*self.objects[i].position[2]
                m=m+self.objects[i].mass
            self.mass_center=[xc/m,yc/m,zc/m]
            self.mass=m
            
    def forces(self,objectI,theta):
        if self.n_objects>0:
            width=self.a-self.b
            v = self.mass_center - objectI.position
            d = sum((v)**2) ** 0.5
            #somehow sometimes it applies force to itself if leave check d!=0
            if d>1e-5:
                t=width/d
                if t < theta or self.children[0]==None:
                    m=self.mass
                    a = LOM_Force(1, m, d)*v/d
                    return a
                else:
                    a1=self.children[0].forces(objectI,theta)
                    a2=self.children[1].forces(objectI,theta)
                    a3=self.children[2].forces(objectI,theta)
                    a4=self.children[3].forces(objectI,theta)
                    a5=self.children[4].forces(objectI,theta)
                    a6=self.children[5].forces(objectI,theta)
                    a7=self.children[6].forces(objectI,theta)
                    a8=self.children[7].forces(objectI,theta)
                    return a1+a2+a3+a4+a5+a6+a7+a8
            else:
                return np.array([0,0,0])
        else:
            return np.array([0,0,0])