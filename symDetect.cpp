#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <list>
#include <utility>
#include <string>
#include <cstdio>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <regex>
#include <tuple>



/// to do
// fix sampling start there
// figure out rho stuff 

using namespace std;
class Point
{
    public:
        int index;
        int x,y;
        // index of each edge
        int edge1=-1,edge2=-1;
        double curvature;
        double tangent;
        // corisponding weight/length fo reach edge
        Point(int i,int xIn, int yIn)
        {   
            index = i;
            x =xIn,y=yIn;
        }
        double sqr_dist(Point p)
        {
            return sqrt(pow(x-p.x,2) + pow(y-p.y,2));
        }
        void addEdge(int p)
        {
            if(edge1>0)
                edge2 = p;
            else
                edge1 = p;
        }
};
struct symPoint
{
    double rho,theta;
    int p1,p2;
    bool near = true;
};
tuple<vector<Point>,vector<int>> buildGraph(string filename)
{
    vector<Point> vertexs,orderedVertex;
    vector<int> samplePoints;
    ifstream obj(filename);
    if(!obj)
    {
        cout<<"Invalid file"<<endl;
        return make_tuple(vertexs, samplePoints);
    }
    // reading from obj file just to keep it simple
    int i = 0;
    string line;
    getline(obj,line);
    // set vertex values
    while(line[0] == 'v')
    {
        int x,y;
        sscanf(line.c_str(),"v %i %i 0",&x,&y);
        vertexs.push_back(Point(i,x,y));
        i+= 1;
        getline(obj,line);
    }
    //sets edge values
    while(!obj.eof())
    {
        int p1,p2;
        sscanf(line.c_str(),"l %i %i",&p1,&p2);
        vertexs[p1-1].addEdge(p2-1);
        vertexs[p2-1].addEdge(p1-1);
        getline(obj,line);
    }
    obj.close();

    // ordering graph
    orderedVertex.push_back(vertexs[0]);
    int prev = 0;
    int c =vertexs[0].edge1;
    int s = 1;
    while(c != -1)
    {
        //vertexs[c].index = s;
        //orderedVertex.push_back(vertexs[c]);
        // butterfly 65
        // star 20
        if(s%((int) vertexs.size()/40) == 0)
        {
            samplePoints.push_back(vertexs[c].index);
        }
        s+=1;
        if(vertexs[c].edge1 != prev)
        {
            prev = c;
            c = vertexs[c].edge1;
        }
        else
        {
            prev = c;
            c = vertexs[c].edge2;
        }
    }

  

    return make_tuple(vertexs, samplePoints);
}

//currently unsued led to weird slope values
double lagrange(vector<int>x,vector<int>y, double z)
{
    double val = 0;
    int size = x.size();
    for(int i= 0; i < size;i++) 
    {
        double term = 1;
        for(int j= 0; j < size;j++)
        {
            if(j!=i)
                term *= term*((z-x[j])/(x[i]-x[j]));
        }
        val += y[i]* term;
    }
    return val;
}
double distanceCalc(Point x, Point y)
{
    return sqrt((pow(x.x - y.x,2) + pow(x.y - y.y,2)));
}
// calculates tangent at point
// uses legrange interpolation and numerical differentiaion using 2 nearest points which have a slope non 0 from given starting point 
void tangentAtPoint(vector<int> sP, vector<Point> &v)
{
    // for dif precision

  //  double h = 0.01;
    // gotta rework    
    for(auto t:sP)
    {
        Point s = v[t];
        int s1 = s.edge1,s2 = s.edge2;
        int prev1 = t, prev2 = t;
        
        while(v[s1].x == s.x || v[s1].y == s.y)
        {
            
            if(v[s1].edge1 == prev1)
            {
                prev1 = s1; 
                s1 = v[s1].edge2;
            }
            else
            {
                prev1 = s1; 
                s1 = v[s1].edge1;
            }
            
        }
        while(v[s2].x == s.x  || v[s2].y == s.y )
        {
            if(v[s2].edge1 == prev2)
            {
                prev2 = s2; 
                s2 = v[s2].edge2;
            }
            else
            {
                prev2 = s2; 
                s2 = v[s2].edge1;
            }
        }

        //changing to menger curve instead of tang for paring
        double a = distanceCalc(s,v[s1]), b = distanceCalc(s,v[s2]), c = distanceCalc(v[s1],v[s2]);
        double p = (a+b+c)/2;
        //hero's formula
        double Area = sqrt(p*(p-a)*(p-b)*(p-c));
        double curvature = (4*Area)/(a*b*c);
        if(curvature < 0.01)
            v[t].curvature = 0;
        else
            v[t].curvature = curvature;
        
        double tanSlope;
        if(v[s1].x >= v[s2].x)              
            v[t].tangent = (double)(v[s1].y-v[s2].y)/( (abs(s.x - v[s2].x)+abs(v[s1].x - s.x))/2);
        else
            v[t].tangent = (double)(v[s2].y-v[s1].y)/( (abs(s.x - v[s2].x)+abs(v[s1].x - s.x))/2);
        // add back tangent as well
        /*
        cout<<"Curvature: "<<s.curvature<<endl;
        
        cout<<endl;
        if(v[s1].x >= v[s2].x)
        {
            cout<<(double)(v[s1].y-v[s2].y)/( abs(v[s1].x - s.x)+abs(v[s2].x - s.x))<<endl;

        } 
        else
            cout<<(double)(v[s2].y-v[s1].y)/( abs(v[s1].x - s.x)+abs(v[s2].x - s.x))<<endl;
        cout<<endl;    */    
    }
}
void writeSample(vector<int> sP, vector<Point> &v,vector<symPoint> m,vector<symPoint> Prem)
{
    ofstream ObjFile,tF,pF;
    tF.open("rhoPoints.txt");
    pF.open("PrerhoPoints.txt");
	ObjFile.open("sample.obj");
    for (auto i : sP)
    {
        ObjFile << "v" << " " << v[i].x << " " << v[i].y << " " << 0 << endl;

    }
    
    for (auto i : m)
    {
        if(i.near) 
            tF << i.rho << "," << i.theta <<endl;

    }
    for (auto i : Prem)
    {

        pF << i.rho << "," << i.theta <<endl;

    }
    ObjFile << "\n";
	ObjFile.close();
    tF.close();
    pF.close();
}
bool notfind(vector<int> d, int i)
{
    for(auto j: d)
    {
        if(j== i)
            return false;
    }
    return true;
}
vector<symPoint> pairing(vector<int> sP, vector<Point> &v)
{
    // rho,theta pair
    vector<symPoint> transForm;
    vector<pair<double,double>> test;
    vector<int> done;
    ofstream tf;
    tf.open("tanPoints.txt");
    
    for(auto i : sP)
    {
        
        for(auto j : sP)
        {
            // main change is removing the notfind check
            //added in tangent check along with curvature
            if(i != j && abs(v[i].curvature)-abs(v[j].curvature) <0.02 && abs(v[i].tangent +v[j].tangent) <= 0.2
                && notfind(done,j))//abs(v[i].curvature +v[j].curvature) <= 0.2 && notfind(done,j)) 
            {
                
                double mX = (double)(v[i].x+v[j].x)/2,mY = (double)(v[i].y+v[j].y)/2;
                tf <<v[i].x<<","<<v[i].y<<","<<mX<<","<<mY<<","<<v[j].x<<","<<v[j].y<<endl;
                
                //printf("%i: %f, %i: %f\n",i,v[i].curvature,j,v[j].curvature);
                // if slope is not undefined
                if(v[i].x != v[j].x)
                {
                    double tanSlope;
                    if(v[i].x >= v[j].x)              
                        tanSlope = (double)(v[i].y-v[j].y)/( abs(v[i].x - v[j].x));
                    else
                        tanSlope = (double)(v[j].y-v[i].y)/( abs(v[i].x - v[j].x));
                        
                    if(abs(tanSlope) < 0.001)
                        tanSlope = 0;
                    
                    
                    // if slope of normal exists / tanSlope != 0
                    if(tanSlope != 0)
                    {
                        double norSlope = -1/tanSlope;
                        double b = mY -norSlope*mX;
                        double xInt = -b/norSlope;
                        pair<double,double> iV,jV;
                        // you have a pic explaining this
                        iV.first = -1, jV.first = -1;
                        iV.second = 0, jV.second = -norSlope;
                        double angle =  acos((iV.first*jV.first + iV.second*jV.second)/(sqrt(pow(iV.first,2)+pow(iV.second,2))*sqrt(pow(jV.first,2)+pow(jV.second,2))));
                        double theta = 1.570796326790 - angle;
                        if(theta < 0.001)
                            theta = 0;
                        double rho = abs(xInt)*cos(theta);
                        double pi = 3.14159265359;
                        //theta = theta*(180/pi);
                        if(norSlope > 0)
                        {
                            if(xInt >=0)
                            {
                                //transForm.push_back(make_pair(rho,2*(3.14159265359)-theta));
                                //transForm.push_back(make_pair(rho,-theta));
                                symPoint temp;
                                temp.rho = rho;
                                temp.theta = -theta;
                                temp.p1 = i;
                                temp.p2 = j;
                                transForm.push_back(temp);
                            }
                            else
                            {
                                //transForm.push_back(make_pair(rho,(3.14159265359)-theta));
                                //transForm.push_back(make_pair(rho,180-theta));
                                symPoint temp;
                                temp.rho = rho;
                                temp.theta = (3.14159265359)-theta;
                                temp.p1 = i;
                                temp.p2 = j;
                                transForm.push_back(temp);
                                
                            }
                        }
                        else
                        {
                            //transForm.push_back(make_pair(rho,theta));
                            symPoint temp;
                            temp.rho = rho;
                            temp.theta = theta;
                            temp.p1 = i;
                            temp.p2 = j;
                            transForm.push_back(temp);
                        }
                        
                    }
                    else // tanSlope = 0 so norSlope is undefined
                    {
                        //transForm.push_back(make_pair(abs(mX),0));
                        symPoint temp;
                        temp.rho = abs(mX);
                        temp.theta = 0;
                        temp.p1 = i;
                        temp.p2 = j;
                        transForm.push_back(temp);
                        
                    }
                    
                }
                else    // if tanslope is undefined norSlope is 0 & angle of normal after transform is 90 as it is horizontal
                {
                    //transForm.push_back(make_pair(abs(mY),1.570796326790));
                    symPoint temp;
                    temp.rho = abs(mY);
                    temp.theta = 1.570796326790;
                    temp.p1 = i;
                    temp.p2 = j;
                    transForm.push_back(temp);
                     
                }
                
            }
        }
        done.push_back(i);
        
   }
    tf.close();

    return transForm;
}
pair<pair<double,double>,pair<double,double>> normalize(vector<symPoint> &p)
{
    double tMin = p[0].theta, tMax = p[0].theta;
    double rMin = p[0].rho,rMax = p[0].rho;
    for(int i =1; i<p.size();i++)
    {
        if(p[i].rho < rMin)
            rMin = p[i].rho;
        else if(p[i].rho > rMax)
            rMax = p[i].rho;
        if(p[i].theta < tMin)
            tMin = p[i].theta;
        else if(p[i].theta > tMax)
            tMax = p[i].theta;
    }
    for(int i =0; i<p.size();i++)
    {
        p[i].rho = (p[i].rho - rMin)/(rMax-rMin);
        p[i].theta = (p[i].theta - tMin)/(tMax-tMin);
    }
    return make_pair(make_pair(rMin,rMax),make_pair(tMin,tMax));

}
double symDistance(symPoint x,symPoint y)
{
    return sqrt((pow(x.rho - y.rho,2) + pow(x.theta - y.theta,2)));
}

double epanKernal(double u)
{
    return (0.75)*(1-pow(u,2));
}

void meanshift(vector<symPoint> &p)
{
    int iterations = 5;
    /// 0.15 is ight for dist
    double dist = 0.15;
    double h =0.15;
    for(int it =0; it < iterations;it++)
    {
        for(int i =0; i<p.size();i++)
        {
            double nX =0,nY=0,d= 0;
            for(int j =0; j<p.size();j++)
            {
                double distance = symDistance(p[i],p[j]);
                //cout<<distance<<endl;
                if(i!= j && distance < dist)
                {
                    //cout<<distance<<endl;
                    double w = epanKernal(distance/h);
                    nX += w*p[j].rho;
                    nY += w*p[j].theta;
                    d += w;
                }
            }
           if(d != 0.0)
           {
                //cout<<i<<endl;
                symPoint temp;
                temp.rho = nX/d;
                temp.theta = nY/d;
                temp.p1 = p[i].p1;
                temp.p2 = p[i].p2;
                p[i] = temp;
           }
           else
           {
                 
                p[i].near = false;
           }
           
        }
    }
}

int main(int argc, char **argv)
{       
        int d = strlen(argv[1]);
        char* filename = argv[1];
               
        //vector<Point> 
        auto[vertexs,samplePoints] = buildGraph(filename);
        if(vertexs.size() == 0)
            return 0;
  
   
        tangentAtPoint(samplePoints,vertexs);
        vector<symPoint> p = pairing(samplePoints,vertexs);
        cout<<p.size()<<endl;
        pair<pair<double,double>,pair<double,double>> norms = normalize(p);
        vector<symPoint> preP = p;
        meanshift(p);
        cout<<p.size()<<endl;

        writeSample(samplePoints,vertexs,p,preP);        
       
       
        return 0;
}

