//#include <igl/cotmatrix.h>
//#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/mat.hpp>
#include <list>
#include <utility>
#include <string>
#include <cstdio>
#include <cmath>
//#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <regex>
#include <tuple>



using namespace std;

class Point
{
    public:
        int index;
        int x,y;
        // index of each edge
        int edge1=-1,edge2=-1;
        int sampEdge1 = -1, sampEdge2 = -1;
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
    
    double xM, yM;
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


    // ordering sample graph
    int prevSamp = 0;
    while(c != -1)
    {
        //vertexs[c].index = s;
        //orderedVertex.push_back(vertexs[c]);
        // butterfly 65
        // star 20
        if(s%((int) vertexs.size()/40) == 0)
        {
            vertexs[c].sampEdge1 = prevSamp;
            vertexs[prevSamp].sampEdge2 = c;
            prevSamp = c;
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
void writeSample(vector<int> sP, vector<Point> &v,vector<pair<double,double>> m,vector<symPoint> Prem)
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
    
        tF << i.first << "," << i.second <<endl;

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
                                temp.xM = mX;
                                temp.yM = mY;
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
                                temp.xM = mX;
                                temp.yM = mY;
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
                            temp.xM = mX;
                            temp.yM = mY;
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
                        temp.xM = mX;
                        temp.yM = mY;
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
                    temp.xM = mX;
                    temp.yM = mY;
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

vector<pair<double,double>> meanshift(vector<symPoint> &p)
{
    int iterations = 5;
    /// 0.15 is ight for dist
    double dist = 0.15;
    // less but tight
    //double dist = 0.045;
    double h =0.15;
    //performs mean shift
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
                //temp.p1 = p[i].p1;
                //temp.p2 = p[i].p2;
                p[i] = temp;
           }
           else
           {
                 
                p[i].near = false;
           }
           
        }
    }
    // identify clusters
    // section off cluster points, then take avaerage at each cluster to determine "peak"
    vector<vector<pair<double,double>>> cTemp;
    cTemp.push_back(vector<pair<double,double>> {make_pair(p[0].rho,p[0].theta)});
    for(auto i : p)
    {
        if(i.near)
        {
            bool inC = false;
            for(int j = 0; j<cTemp.size();j++)
            {
                // voted on the same cluster max, ie extremely close after mean shift
                if(sqrt((pow(i.rho - cTemp[j][0].first,2) + pow(i.theta - cTemp[j][0].second,2)))<0.01)
                {
                    cTemp[j].push_back(make_pair(i.rho,i.theta));
                    inC = true;
                    break;
                }
            }
            if(!inC)
            {
                cTemp.push_back(vector<pair<double,double>> {make_pair(i.rho,i.theta)});
            }
        }
    }
    vector<pair<double,double>> clusters;
    for(auto i : cTemp)
    {
        // avoid small but distant clusters, may be unnecassary
        if(i.size() > 5)
        {
            double sumX = 0;
            double sumY = 0;
            for(auto j : i)
            {
                sumX +=j.first;
                sumY += j.second;
            }
            clusters.push_back(make_pair(sumX/i.size(),sumY/i.size()));
        }
    }

    return clusters;
}
vector<vector<symPoint>> sigPointsFinder(int &sig, vector<pair<double,double>>clusters, vector<symPoint> points,vector<Point> verts)
{
    // finds transform points that lie wiethin the kernal of a voted on max point
    //vector<symPoint> sigPoints;
    vector<vector<symPoint>> sigClusters;
    ofstream sigP,sigLine;
 
    int largestCluster = 0;
    int ind = 0;
    sigP.open("sigPoints.txt");
    sigLine.open("sigLine.txt");

    double d = 0.05;    // distance from cluster max
    for(auto c : clusters)
    {
        symPoint t;
        t.rho = c.first;
        t.theta = c.second;
        vector<symPoint> temp;
          // easy spot for gpu stuff to add
        for(auto p : points)
        {
            if(symDistance(t,p) <= d)
            {
                //sigPoints.push_back(p);
                temp.push_back(p);
                //sigP << p.xM <<" "<< p.yM <<endl;
            }
        }
        // keeps track of cluster with largest content of points
        cout<<temp.size()<<endl;
        cout<<largestCluster<<endl;

        cout<<endl;
        if(temp.size() > largestCluster)
        {
            largestCluster = temp.size();
            ind = sigClusters.size();
        }
        sigClusters.push_back(temp);
    }
    
    for(auto p : sigClusters[ind] )
    {
        sigP << p.xM <<" "<< p.yM <<endl;
        sigLine << verts[p.p1].x<<","<< verts[p.p1].y<<"/"<< verts[p.p2].x<<","<< verts[p.p2].y<<endl;
    }
    sigP.close();
    sig = ind;
    sigLine.close();
    return sigClusters;
}

tuple<cv::Mat_<double>,cv::Mat_<double>> centroid(list<Point> A, list<Point> B, cv::Mat_<double> &matA, cv::Mat_<double> &matB)
{
    // calculates the centroid and realignes points of each point set
    // dont forget opencv is inverted
    // add zip iterator later to learn some boost
    list<Point>::iterator it1 = A.begin();
    list<Point>::iterator it2 = B.begin();
    cv::Mat_<double> AP(2,A.size()), BP(2,B.size());
    int i = 0;
    //load points into mats
    for(; it1 != A.end() && it2 != B.end(); it1++,it2++)
    {   
        matA(0,i) = it1->x;
        matA(1,i) = it1->y;
        matB(0,i) = it2->x;
        matB(1,i) = it2->y;
        i++;
    }
    cv::Mat col_meanB, col_meanA;
    cv::reduce(matB,col_meanB, 1, cv::REDUCE_AVG);
    cv::reduce(matA,col_meanA, 1, cv::REDUCE_AVG);
    for(int i = 0; i < matA.cols; i++)
    {
        AP.col(i) = matA.col(i) - col_meanA;
        BP.col(i) = matB.col(i) - col_meanB;
    }
    
    return make_tuple(AP,BP);


}
cv::Mat_ <double> rotation(cv::Mat_<double> APrime, cv::Mat_<double> BPrime)
{
    cv::Mat_<double> BTrans;
    cv::transpose(BPrime, BTrans);
    cv::Mat_<double> H = APrime*BTrans;
    auto[U,S,V] = cv::SVD(H);
    cv::transpose(U, U);
    cv::Mat_<double> VT;
    cv::transpose(V, VT);
    cv::Mat_<double> R = VT*U;
    
    if(cv::determinant(R) < 0)
    {
        cout<<"RDET: "<<cv::determinant(R)<<endl;
        cout<<R<<endl;
        
        V.col(1) *= -1;
        
        //cv::transpose(U, U);
        cv::transpose(V, VT);
        R = VT*U;
    }
    
    cout<<R<<endl;
    return R;

}

cv::Mat_ <double> scaling(cv::Mat_<double> APrime, cv::Mat_<double> BPrime)
{
    cv::Mat_<double> Sp= cv::Mat_<double>::zeros(2,1), D= cv::Mat_<double>::zeros(2,1);
   
    for(int i = 0; i < APrime.cols; i++)
    {
        cv::Mat_<double> Temp;
        cv::transpose(BPrime.col(i),Temp);
        D += Temp*BPrime.col(i);
        cv::transpose(APrime.col(i),Temp);
        Sp += Temp*APrime.col(i);
    }
    cv::Mat_<double> T;
    cv::divide(D,Sp,T);
    //cout<<"in scaling: "<<T<<endl;
    //cout<<"in scaling D: "<<D<<endl;
    //cout<<"in scaling SP: "<<Sp<<endl;
    //cv::pow(T,1/2,T);
    T(0,0) = sqrt(T(0,0));
    T(1,0) = sqrt(T(1,0));
    //cout<<"in scaling sqrtT: "<<T<<endl;
    return T;

}

cv::Mat_ <double> Translation(cv::Mat_<double> B, cv::Mat_<double> s,cv::Mat_<double> R,cv::Mat_<double> A)
{
    
    cv::Mat col_meanB, col_meanA;
    cv::reduce(B,col_meanB, 1, cv::REDUCE_AVG);
    cv::reduce(A,col_meanA, 1, cv::REDUCE_AVG);
    cv::Mat_<double> RotA = R*col_meanA;
    
    cv::multiply(s, RotA,RotA);
    

    return col_meanB-RotA;
    


}

double error(cv::Mat_<double> BPrime, cv::Mat_<double> s,cv::Mat_<double> R,cv::Mat_<double> APrime,cv::Mat_<double> t)
{
    //cv::Mat_<double> err = cv::Mat_<double>::zeros(2,1);
    double err = 0;
    cv::Mat_<double> d = cv::Mat_<double>::zeros(2,1);
    cv::Mat_<double> TestTemp = R*APrime;
    cout<<endl;
    cout<<"ERROR"<<endl;
    cout<<"APrime: "<<APrime<<endl;
    cout<<"BPrime: "<<BPrime<<endl;
    cout<<"R: "<<R<<endl;
    //cout<<TestTemp<<endl;
    // bottom row is switching away
    TestTemp.row(0) =s.row(0)*TestTemp.row(0);
    TestTemp.row(1) =s.row(1)*TestTemp.row(1);
    cout<<"t: "<<t<<endl;
    cout<<TestTemp<<endl;
    TestTemp.row(0) =t.row(0)+TestTemp.row(0);
    TestTemp.row(1) =t.row(1)+TestTemp.row(1);

    //TestTemp += t;
    cout<<"TestTemp: "<<TestTemp<<endl;
    cout<<BPrime<<endl;
    TestTemp = abs(BPrime - TestTemp);
    cout<<TestTemp<<endl;
    cout<<"MEANX: "<<cv::mean(TestTemp)<<endl;
    return cv::mean(TestTemp)(0); 
    for(int i = 0; i < BPrime.cols; i++)
    {
        cv::Mat_<double> Temp= (R*APrime.col(i));
        cv::multiply(s,Temp,Temp);
        d = BPrime.col(i) - (Temp + t);
        
        cv::Mat_<double> TempD;
        cv::transpose(d,TempD);
        //cout<<"s:"<<s<<endl;
        //cout<<"R:"<<R<<endl;
        //err += TempD*d;
        //change error to just avg
        //cout<<Temp+t<<endl;
        //cout<<BPrime.col(i)<<endl;
        //cout<<d<<endl;
        cv::pow(d,2,d);
        err += sqrt(d(0,0)+d(1,0));
        //cout<<sqrt(d(0,0))//+d(1,0))<<endl;
    }  
    err = err/BPrime.cols;
    //cout<<err<<endl;
    //cv::sqrt(err,err);
    
    return err;

}
bool patchGrowth (list<Point> &A, list<Point> &B, vector<int> &visited, vector<Point> verts)
{
    // dont forget opencv is inverted
    // if list size is only 2, skip first point
    // if 2 d gives weird results add z in and set to 0
    cv::Mat_<double> matA(2,A.size()), matB(2,B.size());
    auto[APrime, BPrime] = centroid(A, B, matA, matB);
    // rotation calc
    cv::Mat_<double> R = rotation(APrime,BPrime);    
    cv::Mat_<double> s = scaling(APrime,BPrime);
    cv::Mat_<double> t = Translation(matB,s,R,matA);
    // figure out error boy
    //cv::Mat_<double> 
    double err = error(matB,s,R,matA,t);
    
    //double m = cv::mean(err)(0);
    //cout<<"Error: "<<m<<endl;
    if(err <= 5 )
    {
        Point endA = A.back(), endB = B.back();
        if(visited[endA.sampEdge1] == 0)
        {
            A.push_back(verts[endA.sampEdge1]);
        }
        else if(visited[endA.sampEdge2] == 0)
        {
            A.push_back(verts[endA.sampEdge2]);
        }
        else
        {
            return true;
        }
        

        if(visited[endB.sampEdge1] == 0)
        {
            B.push_back(verts[endB.sampEdge1]);
        }
        else if(visited[endB.sampEdge2] == 0)
        {
            B.push_back(verts[endB.sampEdge2]);
        }
        else
        {
            cout<<"BRET"<<endl;
            
            return true;
        }
        patchGrowth(A,B,visited,verts);
        cout<<"rec"<<endl;
        return true;
    }
    cout<<"poping"<<endl;
    cout<<A.size()<<endl;
    cout<<B.size()<<endl;
    A.pop_back();
    B.pop_back();
    return false;
}

vector<pair<double,double>> verification(vector<vector<symPoint>> cluster,vector<Point> verts,ofstream &sets,int sig)
{
    // *********************dont forget opencv is inverted***************
    // ICP links
    // http://www.sci.utah.edu/~shireen/pdfs/tutorials/Elhabian_ICP09.pdf
    // http://nghiaho.com/?page_id=671

    // this is a bit wasteful
    vector<int> visited(verts.size(),0);
    vector<pair<double,double>> patch;
    int t = 0;
    //for(auto c : cluster)
    //{
        vector<symPoint> c = cluster[sig];
        sets<<"PairS"<<endl;
        cout<<t<<endl;
        for(auto sP : c)
        {
            // check they havent been visited yet
            if(visited[sP.p1] == 0 && visited[sP.p2] == 0)
            {
                Point pI = verts[sP.p1],pJ = verts[sP.p2];
                sets<<"m:"<<pI.x<<","<<pI.y<<","<<pJ.x<<","<<pJ.y<<endl;
                // first circle around the points
                //vector<Point> atemp = {verts[pI.sampEdge1], pI, verts[pI.sampEdge2]};
                //vector<Point> btemp = {verts[pJ.sampEdge1], pJ, verts[pJ.sampEdge2]};
                list<Point> aB = {verts[pI.sampEdge1],pI, verts[pI.sampEdge2]};
                list<Point> bB = {verts[pJ.sampEdge1],pJ, verts[pJ.sampEdge2]};
                list<Point> aF = {verts[pI.sampEdge2],pI, verts[pI.sampEdge1]};
                list<Point> bF = {verts[pJ.sampEdge2],pJ, verts[pJ.sampEdge1]};
                visited[verts[pI.sampEdge1].index] = 1;
                visited[verts[pI.sampEdge2].index] = 1;
                visited[verts[pJ.sampEdge1].index] = 1;
                visited[verts[pJ.sampEdge2].index] = 1;

                // divide and conquor
                // split up and check both directions from points seperatly
                // swapping the bs
                bool B = patchGrowth(aB,bB,visited,verts);
                bool F = patchGrowth(aF,bF,visited,verts);
                cout<<"TEAT: "<<t<<"######################################"<<endl;
                if(B)
                {   

                    //sets<<"PairS"<<endl;
                    cout<<"true"<<endl;
                    sets<<"F"<<endl;
                    for(auto i : aB)
                    {
                        sets<<i.x<<","<<i.y<<endl;
                    }
                    sets<<"PairE"<<endl;
                    sets<<"B"<<endl;
                    for(auto i : bB)
                    {
                        sets<<i.x<<","<<i.y<<endl;
                    }

                    sets<<"PairE"<<endl;
                }
                if(F)
                {
                    sets<<"F"<<endl;
                    //cout<<"true"<<endl;
                    for(auto i : aF)
                    {
                        sets<<i.x<<","<<i.y<<endl;
                    }
                    sets<<"PairE"<<endl;
                    sets<<"B"<<endl;
                    for(auto i : bF)
                    {
                        sets<<i.x<<","<<i.y<<endl;
                    }
                    sets<<"PairE"<<endl;
                }
            }
            
        }
        
    //}
    
    return patch;
}

int main(int argc, char **argv)
{      
        int d = strlen(argv[1]);
        char* filename = argv[1];
               
        //vector<Point> 
        auto[vertexs,samplePoints] = buildGraph(filename);
        if(vertexs.size() == 0)
            return 0;
        int c = samplePoints[0];
     
        tangentAtPoint(samplePoints,vertexs);
        vector<symPoint> p = pairing(samplePoints,vertexs);
        cout<<p.size()<<endl;
        pair<pair<double,double>,pair<double,double>> norms = normalize(p);
        vector<symPoint> preP = p;
        vector<pair<double,double>> clusters = meanshift(p);
        cout<<clusters.size()<<endl;
        int sig;
        vector<vector<symPoint>> sigCluster =  sigPointsFinder(sig,clusters,preP,vertexs);
        ofstream sets;
        sets.open("matching.txt");
        verification(sigCluster,vertexs,sets,sig);
        
        sets.close();
        
        writeSample(samplePoints,vertexs,clusters,preP);        
       
       
        return 0;
}

