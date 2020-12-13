#pragma once
#ifndef WALK_HPP_
#define WALK_HPP_

#include "ufunction.hpp"
//#include "EF2d-base.hpp"
#include <cstdlib>
// Not really test
// Just coorect compile error
typedef Mesh2d Mesh;
//typedef Vertex2 Vertex;
typedef Simplex Triangle;

typedef double R;

int kthrough=0,kfind=0;
inline int BinaryRand()
{
#ifdef RAND_MAX
  const long HalfRandMax = RAND_MAX/2;
  return rand() <HalfRandMax;
#else
  return rand() & 16384; // 2^14 (
#endif
}


int  WalkInTriangle(const Mesh & Th,int it, double *lambda,
		    R u, R v, R & dt)
{
  const Triangle & T(Th[it]);
  const R2 Q[3]={R2(T[0]),R2( T[1]),R2(T[2])};
//  int i0=Th(T[0]);
//  int i1=Th(T[1]);
//  int i2=Th(T[2]);

  R2 P  = lambda[0]*Q[0]  + lambda[1]*Q[1]  + lambda[2]*Q[2];

  //  cout << " " << u << " " << v ;
  R2 PF = P + R2(u,v)*dt;

  //  couleur(15);MoveTo( P); LineTo( PF);
  R l[3];
  l[0] = det(PF  ,Q[1],Q[2]);
  l[1] = det(Q[0],PF  ,Q[2]);
  l[2] = det(Q[0],Q[1],PF  );
  R Det = l[0]+l[1]+l[2];
  l[0] /= Det;
  l[1] /= Det;
  l[2] /= Det;
  const R eps = 1e-5;
  int neg[3],k=0;
  int kk=-1;
  if (l[0]>-eps && l[1]>-eps && l[2]>-eps)
    {
      dt =0;
      lambda[0] = l[0];
      lambda[1] = l[1];
      lambda[2] = l[2];
    }
  else
    {

      if (l[0]<eps && lambda[0] != l[0]) neg[k++]=0;
      if (l[1]<eps && lambda[1] != l[1]) neg[k++]=1;
      if (l[2]<eps && lambda[2] != l[2]) neg[k++]=2;
      R eps1 = T.mes     * 1.e-5;

      if (k==2) // 2
	{
	  // let j be the vertex beetween the 2 edges
	  int j = 3-neg[0]-neg[1];
	  //
	  R S = det(P,PF,Q[j]);

	  if (S>eps1)
	    kk = (j+1)%3;
	  else if (S<-eps1)
	    kk = (j+2)%3;
	  else if (BinaryRand())
	    kk = (j+1)%3;
	  else
	    kk = (j+2)%3;

	}
      else if (k==1)
	kk = neg[0];
      if(kk>=0)
	{
	  R d=lambda[kk]-l[kk];

	  assert(d);
	  R coef =  lambda[kk]/d;
	  R coef1 = 1-coef;
	  dt        = dt*coef1;
	  lambda[0] = lambda[0]*coef1 + coef *l[0];
	  lambda[1] = lambda[1]*coef1 + coef *l[1];
	  lambda[2] = lambda[2]*coef1 + coef *l[2];
	  lambda[kk] =0;
	}
    }
  int jj=0;
  R lmx=lambda[0];
  if (lmx<lambda[1])  jj=1, lmx=lambda[1];
  if (lmx<lambda[2])  jj=2, lmx=lambda[2];
  if(lambda[0]<0) lambda[jj] += lambda[0],lambda[0]=0;
  if(lambda[1]<0) lambda[jj] += lambda[1],lambda[1]=0;
  if(lambda[2]<0) lambda[jj] += lambda[2],lambda[2]=0;
  return kk;
}
////////////////////////////////////////////// IMPORTANT ////////////////////////////////

// this function gets a u and calculates the value of our function on point lambda dansa triangle it
// U needs to be a
R vP1(const Mesh & Th,int it, double *lambda,const R *U, void * /* numerotation */)
{
    return lambda[0]*U[Th(it,0)] + lambda[1]*U[Th(it,1)] + lambda[2]*U[Th(it,2)];
}

typedef R (*funcEvalEF)(const Mesh & Th,int it, double *lambda,const R *U, void *dataforvP1) ;
int  WalkInTriangleP1(const Mesh & Th,int it, double *lambda,
                    const R *U,const  R *V, R & dt,funcEvalEF efEF,void *dataFE)
{
    R u  = efEF(Th,it,lambda,U,dataFE);//
    R v  = efEF(Th,it,lambda,V,dataFE);
  return WalkInTriangle( Th,it,lambda,u,v,dt);
}


int WalkP1(const Mesh & Th,int& it, R *l,
         const R* U,const R* V, R dt,funcEvalEF efEF,void *dataFE)
{

    int k=0;
    int j;
    while ( (j=WalkInTriangleP1(Th,it,l,U,V,dt,efEF,dataFE))>=0)
        {
            int jj  = j;
            assert( l[j] == 0);
            R a= l[(j+1)%3], b= l[(j+2)%3];
            int itt =  Th[it].SimplexAdj(jj,j);
            if(itt==it || itt <0)  return -1;
            it = itt;
            l[j]=0;
            l[(j+1)%3] = b;
            l[(j+2)%3] = a;
            assert(k++<1000);
        }
    return it;
}

const Triangle *  Find(const Mesh &Th,R2 P, R2 & Phat,bool & outside,const Triangle * tstart, int & it)
{
  int j;
  it=0;
  if ( tstart )
    it =  Th(tstart);
  else  {
    assert(0); // bug not implementa
    /*
    const Vertex * v=quadtree->NearestVertexWithNormal(P);
    if (!v) {
      v=quadtree->NearestVertex(P);
      assert(v);
    }

    it=Contening(v);
    */
 }

  //int its=it,iib=-1,;
  int iit=-1;
  R delta=-1;
  R2 Phatt;
  int k=0;
  kfind++;
  while (1)
    {
   // loop:
      kthrough++;
      if (k++>=1000)
	{
	  assert(k++<1000);
	}
      int kk,n=0,nl[3];

      const Triangle & K=Th[it];

      const R2 & A(K[0]), & B(K[1]), & C(K[2]);
      R l[3]={0,0,0};
      R area2= K.mes*2;
      R eps =  -area2*1e-6;
      l[0] = det(P,B,C);
      l[1] = det(A,P,C);
      l[2] = area2-l[0]-l[1];
      if (l[0] < eps) nl[n++]=0;
      if (l[1] < eps) nl[n++]=1;
      if (l[2] < eps) nl[n++]=2;
      if (n==0)
	{
	  outside=false;
	  Phat=R2(l[1]/area2,l[2]/area2);
	  return &K;
	}
      else if (n==1)
	j=nl[0];
      else
	{
	  kk=BinaryRand() ? 1 : 0;
	  j= nl[ kk ];
	}

      int jj  = j;
      int itt =  K.SimplexAdj(jj,j);

      if(itt==it || itt <0)
	{
	  if ( n==2 )
	    {
	      jj=j= nl[ 1-kk ];
	      itt =  K.SimplexAdj(jj,j);
	      if (itt && itt != it)
		{
		  it=itt;
		  continue;
		}
	    }
	  // projection du point sur la frontiere
	  l[nl[0]]=0;
	  if(n==2) l[nl[1]]=0;
	  R ll=l[0]+l[1]+l[2];
	  Phat=R2(l[1]/ll,l[2]/ll);
	  R2 PQ(K(Phat),P);
	  R dd=(PQ,PQ);
	  if (dd>delta && iit>=0)
	    {
	      Phat=Phatt;
	      return Th.t+iit;}

            int j0=(j+1)%3,i0= Th(K[j0]);
	  int j1=(j+2)%3,i1= Th(K[j1]);
	  int ii=-1,jj;
	  if ( l[j0]> ll/2 ) ii=i0,jj=i1;
	  if ( l[j1]> ll/2 ) ii=i1,jj=i0;
               // cout << ii << " " << jj << " it = " << it << " " << delta << " " << dd <<  endl;
	  //  pour la gestion de la frontiere
	  //  on se promene sur la frontiere
	  /*
	  if (ii>0 && iib != ii )
	    for (int p=BoundaryAdjacencesHead[ii];p>=0;p=BoundaryAdjacencesLink[p])
	      { int e=p/2, ie=p%2, je=2-ie;
		// cout << number(bedges[e][0]) << " " << number(bedges[e][1]) << endl;
		if (! bedges[e].in( vertices+jj))
		  {
		    iib = ii;
		    iit=it;
		    delta=dd;
		    Phatt=Phat;
		    it= BoundaryElement(e,ie);
		    // cout << "  ------ " << it << " " << Phatt <<  endl;
		    goto loop;
		  }
	      }
*/


	  outside=true;
  //  Phat=R2(l[0]*A.x+l[1]*B.x+l[2]*C.x,l[0]*A.y+l[1]*B.y+l[2]*C.y);
	  if (dd>delta && iit>=0)
	    {Phat=Phatt;return Th.t+iit;}
	  else
	    return Th.t+it;
	}
      it=itt;
    }

}
#endif
