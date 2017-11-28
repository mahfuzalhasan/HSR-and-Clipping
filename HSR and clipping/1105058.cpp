#include <iostream>
#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#define pi 2*acos(0.0)
#define rad pi/180.0

using namespace std;

double EyeX,EyeY,EyeZ;
double LookX,LookY,LookZ;
double UpX,UpY,UpZ;
double FovY,Aspect,Near,Far;
int Screen_width,Screen_height;
int R,G,B;
int array[20];
double zB[1800][1800];
int red[1800][1800];
int green[1800][1800];
int blue[1800][1800];
int countB=0;



ofstream outPut("stage1.txt");
ofstream outPut2("stage2.txt");
ofstream outPut3("stage3.txt");
ofstream outPut4("p.txt");




struct point
{
public:

    double x,y,z;
};


bool compare(point a,point b)
{
    return a.y > b.y;
}


point pixel[1800][1800];


class Vector3D
{
public:
    double x,y,z;
    Vector3D(double a, double b, double c)
    {
        x = a;
        y=b;
        z=c;
    }
    Vector3D()
    {

    }

    double dot(Vector3D p)
    {
        return x*p.x + y*p.y + z*p.z;

    }


    Vector3D cross_Multiplication(Vector3D p)
    {

        double a = y*p.z - z*p.y;
        double b = z*p.x - x*p.z;
        double c = x*p.y - y*p.x;
        Vector3D cross = Vector3D(a,b,c);
        return cross;
    }



    void normalize()
    {
        double temp = sqrt(x*x+y*y+z*z);
        x = x / temp;
        y = y / temp;
        z = z / temp;
    }

};


class Triangle
{
public:
        point p[3];
        int r,g,b;
};

queue<Triangle>triPush;


class Matrix
{
    public:
    double matrix[4][4];
    Matrix()
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                if(i==j)
                    matrix[i][j] = 1.0;
                else
                    matrix[i][j] = 0.0;
            }
        }
    }


    void scaling(double s)
    {
        for(int i=0;i<4-1;i++)
        {
            for(int j=0;j<4-1;j++)
            {
                matrix[i][j] = matrix[i][j]*s;
            }
        }
    }


    void setMatrix(Triangle t)
    {
        for(int i=0;i<4-1;i++)
        {
            matrix[0][i] = t.p[i].x;
            matrix[1][i] = t.p[i].y;
            matrix[2][i] = t.p[i].z;
        }
        /*for(int i=0;i<4-1;i++)
                matrix[4-1][i] = 1;
        matrix[4-1][4-1] = 0;*/
    }

    void normalize()
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4-1;j++)
            {
                matrix[i][j] = matrix[i][j]/matrix[4-1][j];
            }
        }
    }







    void printS(int stage,int precision)
    {
       // cout<<"1st1"<<endl;
       if(stage==1){
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    //cout<<matrix[j][i]<<" ";
                    outPut<<std::fixed<<std::setprecision(precision)<<matrix[j][i]<<" ";
                }
               // cout<<endl;
                outPut<<endl;
            }
            outPut<<endl;
        }

        if(stage==2)
        {
            for(int i=0;i<4-1;i++)
            {
                for(int j=0;j<4-1;j++)
                {

                    outPut2<<std::fixed<<std::setprecision(precision)<<matrix[j][i]<<" ";
                }

                outPut2<<endl;
            }
            outPut2<<endl;
        }

        if(stage==3)
        {
            for(int i=0;i<4-1;i++)
            {
                for(int j=0;j<4-1;j++)
                {

                    outPut3<<std::fixed<<std::setprecision(precision)<<matrix[j][i]<<" ";
                }

                outPut3<<endl;
            }
            outPut3<<endl;
        }
    }



};


Matrix productM(Matrix A, Matrix B)
{
    Matrix holder;

        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                holder.matrix[i][j]=0;
                for(int k=0;k<4;k++)
                {
                    holder.matrix[i][j] += A.matrix[i][k]*B.matrix[k][j];
                }
            }
        }
        return holder;
}




float area(point p1,point p2, point p3)
{
    float fl = abs((p1.x*(p2.y-p3.y) + p2.x*(p3.y-p1.y) + p3.x*(p1.y-p2.y))/2.0);
    return fl;
}
/*float area(point p1, point p2, point p3)
{
   return abs((p1.x*(p2.y-p3.y) + p2.x*(p3.y-p1.y) + p3.x*(p1.y-p2.y))/2.0);
}*/
bool isInside2(point p1, point p2, point p3, point p)
{
    float A = area(p1,p2,p3);
    float B = area(p,p2,p3);
    float C = area(p1,p,p3);
    float D = area(p1,p2,p);

    return (A == B+C+D);
}


float determine (point p1, point p2, point p3)
{
    float a = p1.x - p3.x;
    float b = p2.y - p3.y;
    float c = p2.x - p3.x;
    float d = p1.y - p3.y;
    return (a) * (b) - (c) * (d);
}

bool isInside (point pt, point v1, point v2, point v3)
{
    bool b1, b2, b3;


    if(determine(pt, v1, v2) < 0.0f)
    {
        b1 = true;
    }
    else
    {
        b1 =false;
    }

    if(determine(pt, v2, v3) < 0.0f)
    {
        b2=true;
    }
    else
    {
        b2 = false;
    }

    if(determine(pt, v3, v1) < 0.0f)
    {
        b3 = true;
    }
    else
    {
        b3 = false;
    }


    return ((b1 == b2) && (b2 == b3));
}

/*bool isInside(point p1, point p2, point p3, point p)
{  

   float A = area (p1,p2,p3);
   float A1 = area (p,p2,p3);
   float A2 = area (p1, p, p3);
   float A3 = area (p1, p2, p);   
   return (A == A1 + A2 + A3);
}*/

Matrix Rotate(double dx, double dy, double dz, double angle)
{

            Matrix n;
            double cosM, sinM;
            double a = dx*dx;
            double b = dy*dy;
            double c = dz*dz;
            double d = dx*dy;
            double e = dy*dz;
            double f = dz*dx;
            double theta = (angle*pi/180);



            n.matrix[0][0] = a*(1 - cos(theta)) + cos(theta);
            n.matrix[1][0] = d*(1 - cos(theta)) + sin(theta)*dz;
            n.matrix[2][0] = f*(1 - cos(theta)) - dy*sin(theta);

            n.matrix[0][1] = d*(1 - cos(theta)) - dz*sin(theta);
            n.matrix[1][1] = b*(1 - cos(theta)) + cos(theta);
            n.matrix[2][1] = e*(1 - cos(theta)) + dx*sin(theta);

            n.matrix[0][2] = f*(1 - cos(theta)) + dy*sin(theta);
            n.matrix[1][2] = e*(1 - cos(theta)) - sin(theta)*dx;
            n.matrix[2][2] = c*(1 - cos(theta)) + cos(theta);

            return n;
            //mainMatrix.product(n);
    //A.TP();
    //return A;

}



void setTriangle(Matrix mt,Triangle t)
{
    for(int i=0;i<4-1;i++)
    {
        t.p[i].x = mt.matrix[0][i];
        t.p[i].y = mt.matrix[1][i];
        t.p[i].z = mt.matrix[2][i];
    }
}




Matrix vTrans,pTrans;


void viewT()
{
    Vector3D l(LookX - EyeX, LookY - EyeY, LookZ - EyeZ);
    Vector3D upVec(UpX, UpY, UpZ);
    l.normalize();
    Vector3D r = l.cross_Multiplication(upVec);
    r.normalize();

    Vector3D u = r.cross_Multiplication(l);

    Matrix T;
    T.matrix[0][4-1] = -EyeX;
    T.matrix[1][4-1] = -EyeY;
    T.matrix[2][4-1] = -EyeZ;


    vTrans.matrix[0][0] = r.x;
    vTrans.matrix[0][1] = r.y;
    vTrans.matrix[0][2] = r.z;

    vTrans.matrix[1][0] = u.x;
    vTrans.matrix[1][1] = u.y;
    vTrans.matrix[1][2] = u.z;

    vTrans.matrix[2][0] = -l.x;
    vTrans.matrix[2][1] = -l.y;
    vTrans.matrix[2][2] = -l.z;
    //vTrans.product(T);
    vTrans = productM(vTrans,T);
}


void grahamScan(Triangle tTri)
{


    //while(sizeofQ--)
    //{
/*        Triangle tTri = triPush.front();
        triPush.pop();
        triPush.push(tTri);*/
        point tempP[3];
        for(int i=0;i<3;i++)
        {
            tempP[i].x = tTri.p[i].x;
            tempP[i].y = tTri.p[i].y;
            tempP[i].z = tTri.p[i].z;
        }
        sort(tempP,tempP+3,compare);
       // cout<<"sorted points"<<endl;
       // cout<<tempP[0].y<<"  "<<tempP[1].y<<"   "<<tempP[2].y<<endl;


        for(int i = 0;  i < Screen_height; i++)
        {
            for(int j = 0; j < Screen_width; j++)
            {
                if(isInside(pixel[i][j],tempP[0],tempP[1],tempP[2]))
                {
                    countB++;
                    double tY = pixel[i][j].y;
                    double z1,z2,x1,x2;
                    if(tempP[1].y > tY)
                    {
                        z1 = tempP[1].z-((tempP[1].z-tempP[2].z)*((tempP[1].y-tY)/(tempP[1].y-tempP[2].y)));
                        z2 = tempP[0].z-((tempP[0].z-tempP[2].z)*((tempP[0].y-tY)/(tempP[0].y-tempP[2].y)));

                        x1 = tempP[1].x-((tempP[1].x-tempP[2].x)*((tempP[1].y-tY)/(tempP[1].y-tempP[2].y)));
                        x2 = tempP[0].x-((tempP[0].x-tempP[2].x)*((tempP[0].y-tY)/(tempP[0].y-tempP[2].y)));

                        double zD = z1 - (z1-z2)*((x1-pixel[i][j].x)/(x1-x2));
                       // cout<<"zD:"<<zD<<endl;
                        if(zD <= zB[i][j])
                        {
                            zB[i][j] = zD;
                            red[i][j] = tTri.r;
                            green[i][j] = tTri.g;
                            blue[i][j] = tTri.b;
                        }

                    }
                    else if(tempP[0].y  >= tY)
                    {
                        z1 = tempP[0].z-((tempP[0].z-tempP[1].z)*((tempP[0].y-tY)/(tempP[0].y-tempP[1].y)));
                        z2 = tempP[0].z-((tempP[0].z-tempP[2].z)*((tempP[0].y-tY)/(tempP[0].y-tempP[2].y)));

                        x1 = tempP[0].x-((tempP[0].x-tempP[1].x)*((tempP[0].y-tY)/(tempP[0].y-tempP[1].y)));
                        x2 = tempP[0].x-((tempP[0].x-tempP[2].x)*((tempP[0].y-tY)/(tempP[0].y-tempP[2].y)));

                        double zD = z1 - (z1-z2)*((x1-pixel[i][j].x)/(x1-x2));
                        //cout<<"zD:"<<zD<<endl;
                        if(zD <= zB[i][j])
                        {
                            zB[i][j] = zD;
                            red[i][j] = tTri.r;
                            green[i][j] = tTri.g;
                            blue[i][j] = tTri.b;
                        }
                    }
                }
            }
        }


    //}
}

stack<Matrix> push;


void print(queue<Triangle>triPush)
{
     int n = triPush.size();
            while(n--)
            {
                Triangle t = triPush.front();
                triPush.pop();
                for(int i=0;i<3;i++)
                {
                    outPut4<<std::fixed<<std::setprecision(7)<<t.p[i].x<<" "<<t.p[i].y<<" "<<t.p[i].z<<endl;
                }
                outPut4<<"one triangle ends"<<endl;
            }
}

int main()
{

    ifstream fInput("scene.txt");
    fInput>>EyeX>>EyeY>>EyeZ;
    fInput>>LookX>>LookY>>LookZ;
    fInput>>UpX>>UpY>>UpZ;
    fInput>>FovY>>Aspect>>Near>>Far;
    fInput>>Screen_width>>Screen_height;
    fInput>>R>>G>>B;

    for(int i = 0; i < Screen_height; i++)
    {
        for(int j = 0; j < Screen_width; j++)
        {
            red[i][j] = R;
            green[i][j] = G;
            blue[i][j] = B;
        }
    }





    for(int i = 0; i < Screen_height; i++)
    {
        for(int j = 0; j < Screen_width; j++)
        {
            pixel[i][j].x = (-1) + (1/(double)Screen_width)*(2*j + 1);
            pixel[i][j].y = (1) - (1/(double)Screen_height)*(2*i + 1);
            pixel[i][j].z = 1;
        }
    }


    viewT();

    double FovX = FovY * Aspect;
    double t = Near * tan(rad*FovY/2);
    double r = Near * tan(rad*FovX/2);
    Matrix ob;
    ob.matrix[0][0] = Near/r;
    ob.matrix[1][1] = Near/t;
    ob.matrix[2][2] = -(Far + Near)/(Far - Near);
    ob.matrix[2][3] = -(2 * Far * Near)/(Far - Near);
    ob.matrix[3][2] = -1;
    ob.matrix[3][3] = 0;

    pTrans = ob;


    for(int i = 0; i < Screen_height; i++)
    {
        for(int j = 0; j < Screen_width; j++)
        {
            zB[i][j] = 1.0;
        }
    }

    Matrix mainMatrix;
    bool val = true;
    while(val)
    {
        string s;
        fInput>>s;
        if(s=="triangle")
        {
            while(!triPush.empty())
            {
                triPush.pop();
            }
            point p[3];
            double x,y,z;
            int r,g,b;
             Matrix triangle;
            Triangle tri;
            for(int i=0;i<3;i++)
            {
                fInput>>x>>y>>z;
                triangle.matrix[0][i] = x;
                triangle.matrix[1][i] = y;
                triangle.matrix[2][i] = z;

               /* tri.p[i].x=x;
                tri.p[i].y=y;
                tri.p[i].z=z;*/

            }

            fInput>>r>>g>>b;
            tri.r=r;
            tri.g=g;
            tri.b=b;




            for(int i=0;i<4-1;i++)
                triangle.matrix[4-1][i] = 1;
            triangle.matrix[4-1][4-1] = 0;
            //triangle.printCurrentMatrix();
            Matrix temp = mainMatrix;
            //transFormation.printCurrentMatrix();
            //temp.printCurrentMatrix();
            //temp.product(triangle);
            temp = productM(temp,triangle);
            //temp.TP();
            temp.printS(1,3);
            //cout<<endl;


            Matrix tV = vTrans;
            //tV.product(temp);
            tV = productM(tV,temp);
            tV.printS(2,3);


            for(int i=0;i<3;i++)
            {
                tri.p[i].x = tV.matrix[0][i];
                tri.p[i].y = tV.matrix[1][i];
                tri.p[i].z = tV.matrix[2][i];
            }

// CLIPPING STARTS FROM HERE


            triPush.push(tri);
            int s = triPush.size();
            int countZ;
            int countNear;
            int countFar;
            //Vector3D Nclip;
            //Vector3D Nclip;
            //for near plane clip
           while(s--)
            {
                Triangle t = triPush.front();
                triPush.pop();
                countZ = 0;
                for(int i=0;i<3;i++)
                {
                    //cout<<"t.z"<<t.p[i].z<<endl;
                    if( t.p[i].z > -Near)
                    {
                        countZ++;
                        if(countZ==1)
                        {
                            swap(t.p[0],t.p[i]);
                        }
                        if(countZ==2)
                        {
                            swap(t.p[1],t.p[i]);
                        }
                    }
                }
                if(countZ==0)
                {
                    triPush.push(t);
                }
                else if(countZ==1)//one point outside
                {
                    //two Vector3D with one outside point and two inside points
                    Vector3D Nclip1(t.p[1].x-t.p[0].x, t.p[1].y-t.p[0].y, t.p[1].z-t.p[0].z);
                    Vector3D Nclip2(t.p[2].x-t.p[0].x, t.p[2].y-t.p[0].y, t.p[2].z-t.p[0].z);

                    //Now finding coordinates with straight line and near plane

                    //first straight line(point p0 and p1) and near plane co-ordinates

                    Nclip1.normalize();
                    double tt = (-Near-t.p[0].z)/Nclip1.z;
                    point touch1;
                    touch1.x = t.p[0].x+tt*Nclip1.x;
                    touch1.y = t.p[0].y+tt*Nclip1.y;
                    touch1.z = -Near;

                    //secondly straight line(point p2 and p0) and near plane co-ordinates
                    Nclip2.normalize();
                    tt = (-Near-t.p[0].z)/Nclip2.z;
                    point touch2;
                    touch2.x = t.p[0].x+tt*Nclip2.x;
                    touch2.y = t.p[0].y+tt*Nclip2.y;
                    touch2.z = -Near;

                    //two touch point has been found.

                    //now new 1st triangle
                    Triangle t2;
                    t2.p[0].x = touch1.x;
                    t2.p[0].y = touch1.y;
                    t2.p[0].z = touch1.z;

                    t2.p[1].x = t.p[1].x;
                    t2.p[1].y = t.p[1].y;
                    t2.p[1].z = t.p[1].z;

                    t2.p[2].x = t.p[2].x;
                    t2.p[2].y = t.p[2].y;
                    t2.p[2].z = t.p[2].z;

                    t2.r = t.r;
                    t2.g = t.g;
                    t2.b = t.b;

                    triPush.push(t2);

                    //now new 2nd triangle
                    Triangle t3;
                    t3.p[0].x = touch1.x;
                    t3.p[0].y = touch1.y;
                    t3.p[0].z = touch1.z;

                    t3.p[1].x = touch2.x;
                    t3.p[1].y = touch2.y;
                    t3.p[1].z = touch2.z;

                    t3.p[2].x = t.p[2].x;
                    t3.p[2].y = t.p[2].y;
                    t3.p[2].z = t.p[2].z;

                    t3.r = t.r;
                    t3.g = t.g;
                    t3.b = t.b;

                    triPush.push(t3);


                }
                else if(countZ==2)
                {
                    //two point outside -near p[0] && p[1]
                    Vector3D Nclip1(t.p[2].x-t.p[0].x, t.p[2].y-t.p[0].y, t.p[2].z-t.p[0].z);
                    Vector3D Nclip2(t.p[2].x-t.p[1].x, t.p[2].y-t.p[1].y, t.p[2].z-t.p[1].z);


                    //Now finding coordinates with straight line and near plane

                    //first straight line(point p0 and p1) and near plane co-ordinates

                    Nclip1.normalize();
                    double tt = (-Near-t.p[0].z)/Nclip1.z;
                    point touch1;
                    touch1.x = t.p[0].x+tt*Nclip1.x;
                    touch1.y = t.p[0].y+tt*Nclip1.y;
                    touch1.z = -Near;

                    //secondly straight line(point p2 and p0) and near plane co-ordinates
                    Nclip2.normalize();
                    tt = (-Near-t.p[1].z)/Nclip2.z;
                    point touch2;
                    touch2.x = t.p[1].x + tt*Nclip2.x;
                    touch2.y = t.p[1].y + tt*Nclip2.y;
                    touch2.z = -Near;


                    //two touch point has been found.

                    //now new 1st triangle
                    Triangle t3;
                    t3.p[0].x = touch1.x;
                    t3.p[0].y = touch1.y;
                    t3.p[0].z = touch1.z;

                    t3.p[1].x = touch2.x;
                    t3.p[1].y = touch2.y;
                    t3.p[1].z = touch2.z;

                    t3.p[2].x = t.p[2].x;
                    t3.p[2].y = t.p[2].y;
                    t3.p[2].z = t.p[2].z;

                    t3.r = t.r;
                    t3.g = t.g;
                    t3.b = t.b;

                    triPush.push(t3);
                }

            }

            //for far plane

            s = triPush.size();
           // cout<<"s: "<<s<<endl;

            while(s--)
            {
                Triangle t = triPush.front();
                triPush.pop();
                countZ = 0;
                for(int i=0;i<3;i++)
                {
                    if( t.p[i].z < -Far)
                    {
                        countZ++;
                        if(countZ==1)
                        {
                            swap(t.p[0],t.p[i]);
                        }
                        if(countZ==2)
                        {
                            swap(t.p[1],t.p[i]);
                        }
                    }
                }
                if(countZ==0)
                {
                   // cout<<"nothing in far"<<endl;
                    triPush.push(t);
                }
                else if(countZ==1)//one point outside
                {
                    //two Vector3D with one outside point and two inside points
                    Vector3D Nclip1(t.p[1].x-t.p[0].x, t.p[1].y-t.p[0].y, t.p[1].z-t.p[0].z);
                    Vector3D Nclip2(t.p[2].x-t.p[0].x, t.p[2].y-t.p[0].y, t.p[2].z-t.p[0].z);

                    //Now finding coordinates with straight line and near plane

                    //first straight line(point p0 and p1) and near plane co-ordinates

                    Nclip1.normalize();
                    double tt = (-Far-t.p[0].z)/Nclip1.z;
                    point touch1;
                    touch1.x = t.p[0].x+tt*Nclip1.x;
                    touch1.y = t.p[0].y+tt*Nclip1.y;
                    touch1.z = -Far;

                    //secondly straight line(point p2 and p0) and near plane co-ordinates
                    Nclip2.normalize();
                    tt = (-Far-t.p[0].z)/Nclip2.z;
                    point touch2;
                    touch2.x = t.p[0].x+tt*Nclip2.x;
                    touch2.y = t.p[0].y+tt*Nclip2.y;
                    touch2.z = -Far;

                    //two touch point has been found.

                    //now new 1st triangle
                    Triangle t2;
                    t2.p[0].x = touch1.x;
                    t2.p[0].y = touch1.y;
                    t2.p[0].z = touch1.z;

                    t2.p[1].x = t.p[1].x;
                    t2.p[1].y = t.p[1].y;
                    t2.p[1].z = t.p[1].z;

                    t2.p[2].x = t.p[2].x;
                    t2.p[2].y = t.p[2].y;
                    t2.p[2].z = t.p[2].z;

                    t2.r = t.r;
                    t2.g = t.g;
                    t2.b = t.b;

                    triPush.push(t2);

                    //now new 2nd triangle
                    Triangle t3;
                    t3.p[0].x = touch1.x;
                    t3.p[0].y = touch1.y;
                    t3.p[0].z = touch1.z;

                    t3.p[1].x = touch2.x;
                    t3.p[1].y = touch2.y;
                    t3.p[1].z = touch2.z;

                    t3.p[2].x = t.p[2].x;
                    t3.p[2].y = t.p[2].y;
                    t3.p[2].z = t.p[2].z;

                    t3.r = t.r;
                    t3.g = t.g;
                    t3.b = t.b;

                    triPush.push(t3);


                }
                else if(countZ==2)
                {
                    //two point outside -near p[0] && p[1]
                    Vector3D Nclip1(t.p[2].x-t.p[0].x, t.p[2].y-t.p[0].y, t.p[2].z-t.p[0].z);
                    Vector3D Nclip2(t.p[2].x-t.p[1].x, t.p[2].y-t.p[1].y, t.p[2].z-t.p[1].z);


                    //Now finding coordinates with straight line and near plane

                    //first straight line(point p0 and p1) and near plane co-ordinates

                    Nclip1.normalize();
                    double tt = (-Far-t.p[0].z)/Nclip1.z;
                    point touch1;
                    touch1.x = t.p[0].x+tt*Nclip1.x;
                    touch1.y = t.p[0].y+tt*Nclip1.y;
                    touch1.z = -Far;

                    //secondly straight line(point p2 and p0) and near plane co-ordinates

                    Nclip2.normalize();
                    tt = (-Far-t.p[1].z)/Nclip2.z;
                    point touch2;
                    touch2.x = t.p[1].x+tt*Nclip2.x;
                    touch2.y = t.p[1].y+tt*Nclip2.y;
                    touch2.z = -Far;


                    //two touch point has been found.

                    //now new 1st triangle
                    Triangle t3;
                    t3.p[0].x = touch1.x;
                    t3.p[0].y = touch1.y;
                    t3.p[0].z = touch1.z;

                    t3.p[1].x = touch2.x;
                    t3.p[1].y = touch2.y;
                    t3.p[1].z = touch2.z;

                    t3.p[2].x = t.p[2].x;
                    t3.p[2].y = t.p[2].y;
                    t3.p[2].z = t.p[2].z;

                    t3.r = t.r;
                    t3.g = t.g;
                    t3.b = t.b;

                    triPush.push(t3);
                }

            }

            int sizeofQ = triPush.size();
            while(sizeofQ--)
            {
                Triangle tempT = triPush.front();
                triPush.pop();

                tV.setMatrix(tempT);
               // tV.printStage();


                Matrix m = pTrans;
                //m.product(tV);
                m = productM(m,tV);
                m.normalize();
                m.printS(3,7);

               for(int i=0;i<4-1;i++)
                {
                    tempT.p[i].x = m.matrix[0][i];
                    tempT.p[i].y = m.matrix[1][i];
                    tempT.p[i].z = m.matrix[2][i];
                }
               /* for(int i=0;i<3;i++)
                {
                    cout<<"x:"<<tempT.p[i].x<<" "<<"y:"<<tempT.p[i].y<<" "<<"z:"<<tempT.p[i].z<<endl;
                }*/
                grahamScan(tempT);
                triPush.push(tempT);

            }
            //print(triPush);


           /* Matrix m = pTrans;
            m.product(tV);
            m.normalize();
            m.printStage3();*/
            //triangle.printCurrentMatrix();
        }
        else if(s=="translate")
        {
            double x,y,z;
            fInput>>x>>y>>z;
            Matrix translate ;
            translate.matrix[0][4-1] = x;
            translate.matrix[1][4-1] = y;
            translate.matrix[2][4-1] = z;
            //translate.printCurrentMatrix();
          //  mainMatrix.product(translate);
            mainMatrix = productM(mainMatrix,translate);
            //cout<<endl;
            //transFormation.printCurrentMatrix();
            //cout<<endl;
        }
        else if (s=="scale")
        {
            double x,y,z;
            fInput>>x>>y>>z;
            Matrix scale ;
            scale.matrix[0][0] = x;
            scale.matrix[1][1] = y;
            scale.matrix[2][2] = z;
            //scale.printCurrentMatrix();
            //cout<<endl;
          //  mainMatrix.product(scale);
           mainMatrix = productM(mainMatrix,scale);
            //transFormation.printCurrentMatrix();
            //cout<<endl;
        }
        else if(s=="rotate")
        {
            double dx, dy, dz, angle;
            fInput>>angle>>dx>>dy>>dz;
            Matrix p = Rotate(dx,dy,dz,angle);
           // mainMatrix.product(p);
           mainMatrix =  productM(mainMatrix,p);
           /* double angle,x,y,z;
            fInput>>angle>>x>>y>>z;
            Matrix rotate = Rotation(angle,x,y,z);
            //cout<<"rotation matrix"<<endl;
            //rotate.printCurrentMatrix();
            //cout<<endl;
            mainMatrix.product(rotate);*/
            //transFormation.printCurrentMatrix();
            //cout<<endl;
        }
        else if(s=="push")
        {
            push.push(mainMatrix);
        }
        else if(s=="pop")
        {
            Matrix temp = push.top();
            mainMatrix = temp;
            push.pop();
        }
        else if(s == "end"){
          //  cout<<"printing queue"<<endl;
            /*int n = triPush.size();
            while(n--)
            {
                Triangle t = triPush.front();
                triPush.pop();
                for(int i=0;i<3;i++)
                {
                    outPut4<<std::fixed<<std::setprecision(7)<<t.p[i].x<<" "<<t.p[i].y<<" "<<t.p[i].z<<endl;
                }
                outPut4<<"one triangle ends"<<endl;
            }*/
            val=false;
        }
    }

   bitmap_image image(Screen_width,Screen_height);  // Creating an image
    for(int i=0; i<Screen_width; i++)
        for(int j=0; j<Screen_height; j++)
            image.set_pixel(i,j,red[j][i],green[j][i],blue[j][i]);

    image.save_image("out.bmp");
    cout<<countB;
    //fInput.close();
    //outPut.close();
    //outPut2.close();
    //outPut3.close();
    return 0;
}
