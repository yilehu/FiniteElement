#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <string.h>

const double Pi = 3.141592653589793;

struct node_structure
{
	double X;
	double Y;
};

struct elem_structure
{
	int Node[4];
};

/* Discretization */
void Discretization(struct node_structure *Node,double L,double W,int NumL,int NumW,double &dx,double &dy)
{
	int i,j,k;
	dx = L/(double)NumL;
	dy = W/(double)NumW;
	for(j=0;j<=NumW;j++)
	{
		for(i=0;i<=NumL;i++)
		{
			k = i + j*(NumL+1);
			Node[k].X = (double)i * dx;
			Node[k].Y = (double)j * dy;
		}
	}
}

void Connectivity(struct elem_structure *Elem,int NumL,int NumW)
{
	int i,j,k,m,n;
	for(j=0;j<NumW;j++)
	{
		for(i=0;i<NumL;i++)
		{
			k = i +     j*NumL;
			m = i +     j*(NumL+1);
			n = i + (j+1)*(NumL+1);
			Elem[k].Node[0] = m;
			Elem[k].Node[1] = m+1;
			Elem[k].Node[2] = n+1;
			Elem[k].Node[3] = n;
			//printf("ElemID = %2d, i = %3d, j = %3d, k = %3d, l = %3d\n",k,m,m+1,n+1,n);
		}
	}
}

/* Material */
void ElasticMatrixIsotropic(double **D,double E,double v)
{
	double temp1;
	temp1 = E / (1.0 - v * v);
	D[0][0] = temp1;
	D[0][1] = temp1 * v;
	D[0][2] = 0.0;
	D[1][0] = temp1 * v;
	D[1][1] = temp1;
	D[1][2] = 0.0;
	D[2][0] = 0.0;
	D[2][1] = 0.0;
	D[2][2] = temp1 * (1.0 - v) / 2.0;
}

void ComputeABDIsotropic(double **D,double **AA,double **BB,double **DD,int DDim,double t)
{
	int i,j;
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			AA[i][j] = D[i][j] * t;
			BB[i][j] = D[i][j] * 0.0;
			DD[i][j] = D[i][j] * t*t*t/12.0;
		}
	}
}

/* Elemental Stiffness */
void InitializeElemMatrix(double ****ElemBIntPoint,double **JacobianIntPoint,double ***ElemMatrix,int DDim,int KDim,int NumIntPoint,int i)
{
	int j,k,l;
	for(j=0;j<DDim;j++)
	{
		for(k=0;k<KDim;k++)
		{
			for(l=0;l<NumIntPoint*NumIntPoint;l++)
			{
				ElemBIntPoint[i][l][j][k]=0.0;
				JacobianIntPoint[i][l] = 0.0;
			}
			for(l=0;l<KDim;l++)
			{
				ElemMatrix[i][k][l] =0.0;
			}
		}
	}
}

void InitializeBTDMatrix(double **BTD,int DDim,int KDim)
{
	int k,l;
	for(k=0;k<KDim;k++)
	{
		for(l=0;l<DDim;l++)
		{
			BTD[k][l] = 0.0;
		}
	}
}

void BNumIntegral(struct node_structure *Node, struct elem_structure *Elem, int i, double ****B, double **J, double **X, double **Y)
{
	int j,k;
	double kexi[2],eta[2];
	kexi[0] = -1.0/sqrt(3.0);
	kexi[1] = 1.0/sqrt(3.0);
	eta[0] = -1.0/sqrt(3.0);
	eta[1] = 1.0/sqrt(3.0);
	for(j=0;j<2;j++)
	{
		for(k=0;k<2;k++)
		{
			J[i][2*j+k] = (((Node[Elem[i].Node[2]].X-Node[Elem[i].Node[0]].X)*(Node[Elem[i].Node[3]].Y-Node[Elem[i].Node[1]].Y)+(Node[Elem[i].Node[1]].X-Node[Elem[i].Node[3]].X)*(Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[0]].Y))+((Node[Elem[i].Node[2]].X-Node[Elem[i].Node[3]].X)*(Node[Elem[i].Node[0]].Y-Node[Elem[i].Node[1]].Y)+(Node[Elem[i].Node[1]].X-Node[Elem[i].Node[0]].X)*(Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[3]].Y))*kexi[j] + ((Node[Elem[i].Node[2]].X-Node[Elem[i].Node[1]].X)*(Node[Elem[i].Node[3]].Y-Node[Elem[i].Node[0]].Y)+(Node[Elem[i].Node[0]].X-Node[Elem[i].Node[3]].X)*(Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[1]].Y))*eta[k]) / 8.0;
			X[i][2*j+k] = 0.25*(1.0-kexi[j])*(1.0-eta[k])*Node[Elem[i].Node[0]].X + 0.25*(1.0+kexi[j])*(1.0-eta[k])*Node[Elem[i].Node[1]].X + 0.25*(1.0+kexi[j])*(1.0+eta[k])*Node[Elem[i].Node[2]].X + 0.25*(1.0-kexi[j])*(1.0+eta[k])*Node[Elem[i].Node[3]].X;
			Y[i][2*j+k] = 0.25*(1.0-kexi[j])*(1.0-eta[k])*Node[Elem[i].Node[0]].Y + 0.25*(1.0+kexi[j])*(1.0-eta[k])*Node[Elem[i].Node[1]].Y + 0.25*(1.0+kexi[j])*(1.0+eta[k])*Node[Elem[i].Node[2]].Y + 0.25*(1.0-kexi[j])*(1.0+eta[k])*Node[Elem[i].Node[3]].Y;
			B[i][2*j+k][0][0] = ((Node[Elem[i].Node[1]].Y-Node[Elem[i].Node[3]].Y)+(Node[Elem[i].Node[3]].Y-Node[Elem[i].Node[2]].Y)*kexi[j]+(Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[1]].Y)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][2] = ((Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[0]].Y)+(Node[Elem[i].Node[2]].Y-Node[Elem[i].Node[3]].Y)*kexi[j]+(Node[Elem[i].Node[0]].Y-Node[Elem[i].Node[3]].Y)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][4] = ((Node[Elem[i].Node[3]].Y-Node[Elem[i].Node[1]].Y)+(Node[Elem[i].Node[0]].Y-Node[Elem[i].Node[1]].Y)*kexi[j]+(Node[Elem[i].Node[3]].Y-Node[Elem[i].Node[0]].Y)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][6] = ((Node[Elem[i].Node[0]].Y-Node[Elem[i].Node[2]].Y)+(Node[Elem[i].Node[1]].Y-Node[Elem[i].Node[0]].Y)*kexi[j]+(Node[Elem[i].Node[1]].Y-Node[Elem[i].Node[2]].Y)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][1] = ((Node[Elem[i].Node[3]].X-Node[Elem[i].Node[1]].X)+(Node[Elem[i].Node[2]].X-Node[Elem[i].Node[3]].X)*kexi[j]+(Node[Elem[i].Node[1]].X-Node[Elem[i].Node[2]].X)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][3] = ((Node[Elem[i].Node[0]].X-Node[Elem[i].Node[2]].X)+(Node[Elem[i].Node[3]].X-Node[Elem[i].Node[2]].X)*kexi[j]+(Node[Elem[i].Node[3]].X-Node[Elem[i].Node[0]].X)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][5] = ((Node[Elem[i].Node[1]].X-Node[Elem[i].Node[3]].X)+(Node[Elem[i].Node[1]].X-Node[Elem[i].Node[0]].X)*kexi[j]+(Node[Elem[i].Node[0]].X-Node[Elem[i].Node[3]].X)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][7] = ((Node[Elem[i].Node[2]].X-Node[Elem[i].Node[0]].X)+(Node[Elem[i].Node[0]].X-Node[Elem[i].Node[1]].X)*kexi[j]+(Node[Elem[i].Node[2]].X-Node[Elem[i].Node[1]].X)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][2][0] = B[i][2*j+k][1][1];
			B[i][2*j+k][2][1] = B[i][2*j+k][0][0];
			B[i][2*j+k][2][2] = B[i][2*j+k][1][3];
			B[i][2*j+k][2][3] = B[i][2*j+k][0][2];
			B[i][2*j+k][2][4] = B[i][2*j+k][1][5];
			B[i][2*j+k][2][5] = B[i][2*j+k][0][4];
			B[i][2*j+k][2][6] = B[i][2*j+k][1][7];
			B[i][2*j+k][2][7] = B[i][2*j+k][0][6];
		}
	}
}

void ComputeBTD(double **BTD,double ****B,double **D,int DDim,int KDim,int i,int j)
{
	int k,l,m;
	for(k=0;k<KDim;k++)
	{
		for(l=0;l<DDim;l++)
		{
			for(m=0;m<DDim;m++)
			{
				BTD[k][l] += B[i][j][m][k] * D[m][l];
			}
		}
	}
}
			
void ComputeElemMatrix1(double ***ElemMatrix,double **BTD,double ****ElemBIntPoint,double **JacobianIntPoint,double t,int DDim,int KDim,int i,int j)
{
	int k,l,m;
	for(k=0;k<KDim;k++)
	{
		for(l=0;l<KDim;l++)
		{
			for(m=0;m<DDim;m++)
			{
				ElemMatrix[i][k][l] += BTD[k][m] * ElemBIntPoint[i][j][m][l] * t * JacobianIntPoint[i][j];
			}
		}
	}
}

void ComputeAC(double **A,double **C,double L,double W)
{
	A[0][0] = 2.0/L; A[0][1] =   0.0; A[0][2] =   0.0; A[0][3] =   0.0;
	A[1][0] =   0.0; A[1][1] =   0.0; A[1][2] =   0.0; A[1][3] = 2.0/W;
	A[2][0] =   0.0; A[2][1] = 2.0/L; A[2][2] = 2.0/W; A[2][3] =   0.0;

	C[0][0] = 4.0/L/L; C[0][1] =     0.0; C[0][2] =     0.0;
	C[1][0] =     0.0; C[1][1] = 4.0/W/W; C[1][2] =     0.0;
	C[2][0] =     0.0; C[2][1] =     0.0; C[2][2] = 8.0/L/W;
}

void ComputeATDC(FILE *fp3,double **ATDC,double **A,double **DD,double **C,int ADim,int DDim,int CDim)
{
	int i,j,k;
	double **ATD;
	ATD = (double**)malloc(ADim*sizeof(double*));
	for(i=0;i<ADim;i++)
		{*(ATD+i) = (double*)malloc(DDim*sizeof(double));}

	for(i=0;i<ADim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			ATD[i][j] = 0.0;
			for(k=0;k<DDim;k++)
			{
				ATD[i][j] += A[k][i] * DD[k][j];
			}
		}
	}

	fprintf(fp3,"Stiffness Blocks ATDC\n\n");
	for(i=0;i<ADim;i++)
	{
		for(j=0;j<CDim;j++)
		{
			ATDC[i][j] = 0.0;
			for(k=0;k<DDim;k++)
			{
				ATDC[i][j] += ATD[i][k] * C[k][j];
			}
			fprintf(fp3," %15.6lf",ATDC[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");
}

void ComputeG(double **G,double kexi,double eta,int DimRow,int DimCol)
{
	int i,j;
	for(i=0;i<DimRow;i++)
	{
		for(j=0;j<DimCol;j++)
		{
			G[i][j] = 0.0;
		}
	}
	G[0][0] = -0.25*(1.0 -  eta); G[0][2] =  0.25*(1.0 -  eta); G[0][4] =  0.25*(1.0 +  eta); G[0][6] = -0.25*(1.0 +  eta);
	G[1][1] = -0.25*(1.0 -  eta); G[1][3] =  0.25*(1.0 -  eta); G[1][5] =  0.25*(1.0 +  eta); G[1][7] = -0.25*(1.0 +  eta);
	G[2][0] = -0.25*(1.0 - kexi); G[2][2] = -0.25*(1.0 + kexi); G[2][4] =  0.25*(1.0 + kexi); G[2][6] =  0.25*(1.0 - kexi);
	G[3][1] = -0.25*(1.0 - kexi); G[3][3] = -0.25*(1.0 + kexi); G[3][5] =  0.25*(1.0 + kexi); G[3][7] =  0.25*(1.0 - kexi);
}

void ComputeQ(double **Q,double kexi,double eta,int DimRow,int DimCol,double L,double W)
{
	Q[0][0] =  0.75 *           kexi *(1.0 - eta); Q[1][0] =  0.75 *  (1.0 - kexi)*         eta ; Q[2][0] =  0.125* ( 4.0 - 3.0*kexi*kexi - 3.0* eta* eta); 
	Q[0][1] =                                   0; Q[1][1] = -0.125*W*(1.0 - kexi)*(1.0-3.0*eta); Q[2][1] = -0.0625* W * (-1.0 - 2.0* eta + 3.0* eta* eta); 
	Q[0][2] =  0.125*L*(1.0-3.0*kexi)*(1.0 - eta); Q[1][2] =                                   0; Q[2][2] =  0.0625* L * (-1.0 - 2.0*kexi + 3.0*kexi*kexi); 
	Q[0][3] = -0.75 *           kexi *(1.0 - eta); Q[1][3] =  0.75 *  (1.0 + kexi)*         eta ; Q[2][3] =  0.125* (-4.0 + 3.0*kexi*kexi + 3.0* eta* eta);
	Q[0][4] =                                   0; Q[1][4] = -0.125*W*(1.0 + kexi)*(1.0-3.0*eta); Q[2][4] =  0.0625* W * (-1.0 - 2.0* eta + 3.0* eta* eta);
	Q[0][5] = -0.125*L*(1.0+3.0*kexi)*(1.0 - eta); Q[1][5] =                                   0; Q[2][5] = -0.0625* L * ( 1.0 - 2.0*kexi - 3.0*kexi*kexi);
	Q[0][6] = -0.75 *           kexi *(1.0 + eta); Q[1][6] = -0.75 *  (1.0 + kexi)*         eta ; Q[2][6] =  0.125* ( 4.0 - 3.0*kexi*kexi - 3.0* eta* eta);
	Q[0][7] =                                   0; Q[1][7] =  0.125*W*(1.0 + kexi)*(1.0+3.0*eta); Q[2][7] = -0.0625* W * ( 1.0 - 2.0* eta - 3.0* eta* eta);
	Q[0][8] = -0.125*L*(1.0+3.0*kexi)*(1.0 + eta); Q[1][8] =                                   0; Q[2][8] =  0.0625* L * ( 1.0 - 2.0*kexi - 3.0*kexi*kexi);
	Q[0][9] =  0.75 *           kexi *(1.0 + eta); Q[1][9] = -0.75 *  (1.0 - kexi)*         eta ; Q[2][9] =  0.125* (-4.0 + 3.0*kexi*kexi + 3.0* eta* eta);
	Q[0][10]=                                   0; Q[1][10]=  0.125*W*(1.0 - kexi)*(1.0+3.0*eta); Q[2][10]=  0.0625* W * ( 1.0 - 2.0* eta - 3.0* eta* eta);
	Q[0][11]=  0.125*L*(1.0-3.0*kexi)*(1.0 + eta); Q[1][11]=                                   0; Q[2][11]= -0.0625* L * (-1.0 - 2.0*kexi + 3.0*kexi*kexi);
}

void ComputeGTDQ(double **ElemMatrix,double **G,double **DDD,double **Q,int DimGRow,int DimGCol,int DimQRow,int DimQCol,int i0,int j0)
{
	int i,j,k;
	double **GTD;
	GTD = (double**)malloc(DimGCol*sizeof(double*));
	for(i=0;i<DimGCol;i++)
		{*(GTD+i) = (double*)malloc(DimQRow*sizeof(double));}

	for(i=0;i<DimGCol;i++)
	{
		for(j=0;j<DimQRow;j++)
		{
			GTD[i][j] = 0.0;
			for(k=0;k<DimGRow;k++)
			{
				GTD[i][j] += G[k][i] * DDD[k][j];
			}
		}
	}

	for(i=0;i<DimGCol;i++)
	{
		for(j=0;j<DimQCol;j++)
		{
			for(k=0;k<DimQRow;k++)
			{
				ElemMatrix[i0+i][j0+j] += GTD[i][k] * Q[k][j];
			}
		}
	}
}

void ComputeElemMatrix2(FILE *fp3,double **ElemMatrix,double **AA,double **BB,double **DD,double **A,double **C,int DDim,int NumNodePerElem,int DoF0,int DoF,double L,double W)
{
	int i,j,DimG[2],DimQ[2];
	double detJ = L*W/4.0;
	DimG[0] = DoF0*DoF0;
	DimG[1] = DoF0*NumNodePerElem;
	DimQ[0] = DDim;
	DimQ[1] = (DoF-DoF0)*NumNodePerElem;

	double **G,**Q;
	G = (double**)malloc(DimG[0]*sizeof(double*));
	for(i=0;i<DimG[0];i++)
		{*(G+i) = (double*)malloc(DimG[1]*sizeof(double));}
	Q = (double**)malloc(DimQ[0]*sizeof(double*));
	for(i=0;i<DimQ[0];i++)
		{*(Q+i) = (double*)malloc(DimQ[1]*sizeof(double));}
	
	double **AAA,**BBB,**DDD;
	AAA = (double**)malloc(DimG[0]*sizeof(double*));
	for(i=0;i<DimG[0];i++)
		{*(AAA+i) = (double*)malloc(DimG[0]*sizeof(double));}
	BBB = (double**)malloc(DimG[0]*sizeof(double*));
	for(i=0;i<DimG[0];i++)
		{*(BBB+i) = (double*)malloc(DimQ[0]*sizeof(double));}
	DDD = (double**)malloc(DimQ[0]*sizeof(double*));
	for(i=0;i<DimQ[0];i++)
		{*(DDD+i) = (double*)malloc(DimQ[0]*sizeof(double));}

	ComputeATDC(fp3,AAA,A,AA,A,DimG[0],DDim,DimG[0]);
	ComputeATDC(fp3,BBB,A,BB,C,DimG[0],DDim,DimQ[0]);
	ComputeATDC(fp3,DDD,C,DD,C,DimQ[0],DDim,DimQ[0]);

	double kexi[2],eta[2],Wi = 1.0,Wj = 1.0;
	kexi[0] = -1.0/sqrt(3.0);kexi[1] = 1.0/sqrt(3.0);
	eta[0] = -1.0/sqrt(3.0);eta[1] = 1.0/sqrt(3.0);
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			ComputeG(G,kexi[i],eta[j],DimG[0],DimG[1]);
			ComputeQ(Q,kexi[i],eta[j],DimQ[0],DimQ[1],L,W);
			ComputeGTDQ(ElemMatrix,G,AAA,G,DimG[0],DimG[1],DimG[0],DimG[1],      0,      0);
			ComputeGTDQ(ElemMatrix,G,BBB,Q,DimG[0],DimG[1],DimQ[0],DimQ[1],      0,DimG[1]);
			ComputeGTDQ(ElemMatrix,Q,DDD,Q,DimQ[0],DimQ[1],DimQ[0],DimQ[1],DimG[1],DimG[1]);
		}
	}

	//double kexi = 0.0,eta = 0.0,Wi = 2.0,Wj = 2.0;
	//ComputeG(G,kexi,eta,DimG[0],DimG[1]);
	//ComputeQ(Q,kexi,eta,DimQ[0],DimQ[1],L,W);
	//ComputeGTDQ(ElemMatrix,G,AAA,G,DimG[0],DimG[1],DimG[0],DimG[1],      0,      0);
	//ComputeGTDQ(ElemMatrix,G,BBB,Q,DimG[0],DimG[1],DimQ[0],DimQ[1],      0,DimG[1]);
	//ComputeGTDQ(ElemMatrix,Q,DDD,Q,DimQ[0],DimQ[1],DimQ[0],DimQ[1],DimG[1],DimG[1]);


	for(i=0;i<DimG[1];i++)
	{
		for(j=0;j<DimQ[1];j++)
		{
			ElemMatrix[i+0][j+DimG[1]] *= -1.0;
			ElemMatrix[j+DimG[1]][i+0] = ElemMatrix[i+0][j+DimG[1]];
		}
	}

	for(i=0;i<DoF*NumNodePerElem;i++)
	{
		for(j=0;j<DoF*NumNodePerElem;j++)
		{
			ElemMatrix[i][j] *= detJ*Wi*Wj;
		}
	}
}

/* Global Stiffness and Equivalent Nodal Force */
void ComputeGlobalMatrix1(double **GlobMatrix,double ***ElemMatrix,struct elem_structure *Elem,int ElemNum,int DoF,int NumNodePerElem)
{
	int i,j,k;
	for(i=0;i<ElemNum;i++)
	{
		for(j=0;j<NumNodePerElem;j++)
		{
			for(k=0;k<NumNodePerElem;k++)
			{
				GlobMatrix[DoF*Elem[i].Node[j]+0][DoF*Elem[i].Node[k]+0] += ElemMatrix[i][DoF*j+0][DoF*k+0];
				GlobMatrix[DoF*Elem[i].Node[j]+0][DoF*Elem[i].Node[k]+1] += ElemMatrix[i][DoF*j+0][DoF*k+1];
				GlobMatrix[DoF*Elem[i].Node[j]+1][DoF*Elem[i].Node[k]+0] += ElemMatrix[i][DoF*j+1][DoF*k+0];
				GlobMatrix[DoF*Elem[i].Node[j]+1][DoF*Elem[i].Node[k]+1] += ElemMatrix[i][DoF*j+1][DoF*k+1];
			}
		}
	}
}

void ComputeGlobalMatrix2(double **GlobMatrix,double **ElemMatrix,struct elem_structure *Elem,int ElemNum,int NodeNum,int DoF0,int DoF,int NumNodePerElem)
{
	int i,j,k;
	int offsetGlobalx,offsetGlobaly,offsetElemx,offsetElemy;
	offsetGlobalx = NodeNum * DoF0;
	offsetGlobaly = NodeNum * DoF0;
	offsetElemx = NumNodePerElem * DoF0;
	offsetElemy = NumNodePerElem * DoF0;

	for(i=0;i<ElemNum;i++)
	{
		for(j=0;j<NumNodePerElem;j++)
		{
			for(k=0;k<NumNodePerElem;k++)
			{
				GlobMatrix[DoF0*Elem[i].Node[j]+0][DoF0*Elem[i].Node[k]+0] += ElemMatrix[DoF0*j+0][DoF0*k+0];
				GlobMatrix[DoF0*Elem[i].Node[j]+0][DoF0*Elem[i].Node[k]+1] += ElemMatrix[DoF0*j+0][DoF0*k+1];
				GlobMatrix[DoF0*Elem[i].Node[j]+1][DoF0*Elem[i].Node[k]+0] += ElemMatrix[DoF0*j+1][DoF0*k+0];
				GlobMatrix[DoF0*Elem[i].Node[j]+1][DoF0*Elem[i].Node[k]+1] += ElemMatrix[DoF0*j+1][DoF0*k+1];

				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+0][offsetElemy + (DoF-DoF0)*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+0][offsetElemy + (DoF-DoF0)*k+1];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+2] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+0][offsetElemy + (DoF-DoF0)*k+2];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+1][offsetElemy + (DoF-DoF0)*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+1][offsetElemy + (DoF-DoF0)*k+1];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+2] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+1][offsetElemy + (DoF-DoF0)*k+2];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+2][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+2][offsetElemy + (DoF-DoF0)*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+2][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+2][offsetElemy + (DoF-DoF0)*k+1];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+2][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+2] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+2][offsetElemy + (DoF-DoF0)*k+2];

				GlobMatrix[DoF0*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+0] += ElemMatrix[DoF0*j+0][offsetElemy + (DoF-DoF0)*k+0];
				GlobMatrix[DoF0*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+1] += ElemMatrix[DoF0*j+0][offsetElemy + (DoF-DoF0)*k+1];
				GlobMatrix[DoF0*Elem[i].Node[j]+0][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+2] += ElemMatrix[DoF0*j+0][offsetElemy + (DoF-DoF0)*k+2];
				GlobMatrix[DoF0*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+0] += ElemMatrix[DoF0*j+1][offsetElemy + (DoF-DoF0)*k+0];
				GlobMatrix[DoF0*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+1] += ElemMatrix[DoF0*j+1][offsetElemy + (DoF-DoF0)*k+1];
				GlobMatrix[DoF0*Elem[i].Node[j]+1][offsetGlobaly + (DoF-DoF0)*Elem[i].Node[k]+2] += ElemMatrix[DoF0*j+1][offsetElemy + (DoF-DoF0)*k+2];

				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+0][DoF0*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+0][DoF0*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+0][DoF0*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+0][DoF0*k+1];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+1][DoF0*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+1][DoF0*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+1][DoF0*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+1][DoF0*k+1];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+2][DoF0*Elem[i].Node[k]+0] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+2][DoF0*k+0];
				GlobMatrix[offsetGlobalx + (DoF-DoF0)*Elem[i].Node[j]+2][DoF0*Elem[i].Node[k]+1] += ElemMatrix[offsetElemx + (DoF-DoF0)*j+2][DoF0*k+1];
			}
		}
	}


}

void InitializeNodalVarible(double *Nodal_Force,double *Nodal_Displacement,int DoF,int NodeNum)
{
	int i;
	for(i=0;i<DoF*NodeNum;i++)
	{
		Nodal_Force[i] = 0.0;
		Nodal_Displacement[i] = 0.0;
	}
}

void ComputeFeq1(struct node_structure *Node,double *Nodal_Force,int DoF,int NumL,int NumW,double W,double t,double Ta,double Tb)
{
	int i,j,k;
	double *T;
	T = (double*)malloc((NumW+1)*sizeof(double));
	for(i=0;i<NumW+1;i++)
	{
		j = (i+1)*(NumL+1) - 1;
		T[i] = Tb + Node[j].Y * (Ta - Tb) / W;
	}
	for(i=0;i<NumW;i++)
	{
		j = (i+1)*(NumL+1) - 1;
		k = (i+2)*(NumL+1) - 1;
		//Nodal_Force[DoF*j+0] += t*(Node[k].coordy-Node[j].coordy) * (2*T[i+0] + T[i+1])/6.0;
		//Nodal_Force[DoF*k+0] += t*(Node[k].coordy-Node[j].coordy) * (T[i+0] + 2*T[i+1])/6.0;
		Nodal_Force[DoF*j+0] += t*(Node[k].Y-Node[j].Y) * (T[i+0] + T[i+1])/4.0;
		Nodal_Force[DoF*k+0] += t*(Node[k].Y-Node[j].Y) * (T[i+1] + T[i+0])/4.0;
	}
}

void ComputeFeq2(struct elem_structure *Elem,double *Nodal_Force,int DoF0,int DoF,int NodeNum,int NumL,int NumW,double dx,double dy,double q0)
{
	int i,j,k,m,n,p,q;
	
	for(j=0;j<NumW;j++)
	{
		for(i=0;i<NumL;i++)
		{
			k = i + j*NumL;
			m = Elem[k].Node[0];
			n = Elem[k].Node[1];
			p = Elem[k].Node[2];
			q = Elem[k].Node[3];
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*m+0] += q0*dx*dy/4.0;
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*m+1] += q0*dx*dy/4.0 * ( dy/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*m+2] += q0*dx*dy/4.0 * (-dx/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*n+0] += q0*dx*dy/4.0;
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*n+1] += q0*dx*dy/4.0 * ( dy/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*n+2] += q0*dx*dy/4.0 * ( dx/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*p+0] += q0*dx*dy/4.0;
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*p+1] += q0*dx*dy/4.0 * (-dy/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*p+2] += q0*dx*dy/4.0 * ( dx/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*q+0] += q0*dx*dy/4.0;
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*q+1] += q0*dx*dy/4.0 * (-dy/6.0);
			Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*q+2] += q0*dx*dy/4.0 * (-dx/6.0);
		}
		
	}
}

/* Boundary Conditions */
void ComputeBC1(int *BCImplicitStatus,int JDim,int DoF,int NumL,int NumW)
{
	int i,j;
	for(i=0;i<JDim;i++)
	{
		BCImplicitStatus[i] = 0;
	}
	for(i=0;i<NumW+1;i++)
	{
		j = i*(NumL+1);
		BCImplicitStatus[DoF*j+0] = 1;
		BCImplicitStatus[DoF*j+1] = 1;
	}
}

void ComputeBC2(int *BCImplicitStatus,int JDim,int DoF0,int DoF,int NodeNum,int NumL,int NumW)
{
	int i,j;
	for(i=0;i<JDim;i++)
	{
		BCImplicitStatus[i] = 0;
	}
	for(i=0;i<NumW+1;i++)
	{
		j = i*(NumL+1);
		BCImplicitStatus[DoF0*j+0] = 1;
		BCImplicitStatus[DoF0*j+1] = 1;
		BCImplicitStatus[DoF0*NodeNum + (DoF-DoF0)*j+0] = 1;

		j = i*(NumL+1) + NumL;
		BCImplicitStatus[DoF0*j+0] = 1;
		BCImplicitStatus[DoF0*j+1] = 1;
		BCImplicitStatus[DoF0*NodeNum + (DoF-DoF0)*j+0] = 1;
	}
	for(i=0;i<NumL+1;i++)
	{
		j = i;
		BCImplicitStatus[DoF0*j+0] = 1;
		BCImplicitStatus[DoF0*j+1] = 1;
		BCImplicitStatus[DoF0*NodeNum + (DoF-DoF0)*j+0] = 1;

		j = NumW*(NumL+1) + i;
		BCImplicitStatus[DoF0*j+0] = 1;
		BCImplicitStatus[DoF0*j+1] = 1;
		BCImplicitStatus[DoF0*NodeNum + (DoF-DoF0)*j+0] = 1;
	}
}

/* Reduced Stiffness and Equivalent Nodal Force */
int ComputeReducedJDim(int *BCImplicitStatus,int JDim)
{
	int i,ReducedJDim;
	for(i=0,ReducedJDim=0;i<JDim;i++)
	{
		if(BCImplicitStatus[i] == 0)
		{
			ReducedJDim++;
		}
	}
	return ReducedJDim;
}

void ComputeReducedJRow(int *BCImplicitStatus,int *ReducedJRow,int JDim)
{
	int i,j;
	for(i=0,j=0;i<JDim;i++)
	{
		if(BCImplicitStatus[i] == 0)
		{
			ReducedJRow[j] = i;
			j++;
		}
	}
}

void ComputeReducedMatrix(double **GlobMatrix,double *Nodal_Force,double **ReducedMatrix,double *ReducedForce,double *ReducedDisp,int *ReducedJRow,int ReducedJDim)
{
	int i,j;
	for(i=0;i<ReducedJDim;i++)
	{
		ReducedForce[i] = Nodal_Force[ReducedJRow[i]];
		ReducedDisp[i] = 0.0;
		for(j=0;j<ReducedJDim;j++)
		{
			ReducedMatrix[i][j] = GlobMatrix[ReducedJRow[i]][ReducedJRow[j]];
		}
	}
}

/* Solver */
void Inverse(double **R,double **Rinv,double *ReducedForce,double *ReducedDisp,int num)
{
	int i,j,k;
	double temp1;
	double **I;
	I = (double**)malloc(num*sizeof(double*));
	for(i=0;i<num;i++)
	{
		*(I+i) = (double*)malloc(num*sizeof(double));
	}
	for(i=0;i<num;i++)
	{
		for(j=0;j<num;j++)
		{
			Rinv[i][j] = R[i][j];
			if(i==j)
			{
				I[i][j] = 1.0;
			}
			else
			{
				I[i][j] = 0.0;
			}
		}
	}
	for(k=0;k<num-1;k++)
	{
		for(i=k+1;i<num;i++)
		{
			temp1 = Rinv[i][k] / Rinv[k][k];
			for(j=0;j<num;j++)
			{
				Rinv[i][j] -= Rinv[k][j] * temp1;
				I[i][j] -= I[k][j] * temp1;
			}
		}
	}
	for(k=num-1;k>0;k--)
	{
		for(i=k-1;i>-1;i--)
		{
			temp1 = Rinv[i][k] / Rinv[k][k];
			for(j=num-1;j>-1;j--)
			{
				Rinv[i][j] -= Rinv[k][j] * temp1;
				I[i][j] -= I[k][j] * temp1;
			}
		}
	}
	for(i=0;i<num;i++)
	{
		temp1 = Rinv[i][i];
		for(j=0;j<num;j++)
		{
			Rinv[i][j] = I[i][j] / temp1;
		}
	}
	for(i=0;i<num;i++)
	{
		for(j=0;j<num;j++)
		{
			ReducedDisp[i] += Rinv[i][j] * ReducedForce[j];
		}
	}
}

void GaussElimination(double **A,double **U,double *b,int ReducedJDim)
{
	int i,j,k;
	double l;
	for(i=0;i<ReducedJDim;i++)
	{
		for(j=0;j<ReducedJDim;j++)
		{
			U[i][j] = A[i][j];
		}
	}
	for(k=0;k<ReducedJDim-1;k++)
	{
		for(i=k+1;i<ReducedJDim;i++)
		{
			l = U[i][k]/U[k][k];
			for(j=k+1;j<ReducedJDim;j++)
			{
				U[i][j] = U[i][j] - l*U[k][j];
			}
			b[i] = b[i] - l*b[k];
		}
	}
}

void Uxy(double **U,double *x,double *y,int ReducedJDim)
{
	int i,j,k;
	x[ReducedJDim-1] =  y[ReducedJDim-1]/U[ReducedJDim-1][ReducedJDim-1];
	for(k=0;k<ReducedJDim-1;k++)
	{
		i = ReducedJDim-2-k;
		x[i] = y[i];
		for(j=i+1;j<ReducedJDim;j++)
		{
			x[i] -= U[i][j]*x[j];
		}
		x[i] /= U[i][i];
	}
}

void LUDecomposition(double **A,double **LU,int ReducedJDim)
{
	int i,j,k;
	for(i=0;i<ReducedJDim;i++)
	{
		for(j=0;j<ReducedJDim;j++)
		{
			LU[i][j] = A[i][j];
		}
	}
	for(k=0;k<ReducedJDim-1;k++)
	{
		for(j=k;j<ReducedJDim;j++)
		{
			for(i=0;i<k;i++)
			{
				LU[k][j] -= LU[k][i]*LU[i][j];
			}
		}
		for(i=k+1;i<ReducedJDim;i++)
		{
			for(j=0;j<k;j++)
			{
				LU[i][k] -= LU[i][j]*LU[j][k];
			}
			LU[i][k] /= LU[k][k];
		}
	}
	for(i=0;i<ReducedJDim-1;i++)
	{
		LU[ReducedJDim-1][ReducedJDim-1] -= LU[ReducedJDim-1][i]*LU[i][ReducedJDim-1];
	}
}

void Lyb(double **L,double *y,double *b,int ReducedJDim)
{
	int i,j;
	y[0] = b[0];
	for(i=1;i<ReducedJDim;i++)
	{
		y[i] = b[i];
		for(j=0;j<i;j++)
		{
			y[i] -= L[i][j]*y[j];
		}
	}
}

void LLDecomposition(double **A,double **LL,int ReducedJDim)
{
	int i,j,k;
	double sum;
	for(j=0;j<ReducedJDim-1;j++)
	{
		sum = 0.0;
		for(k=0;k<j;k++)
		{
			sum += LL[j][k]*LL[j][k];
		}
		LL[j][j] = sqrt(A[j][j] - sum);
		for(i=j+1;i<ReducedJDim;i++)
		{
			sum = 0.0;
			for(k=0;k<j;k++)
			{
				sum += LL[i][k]*LL[j][k];
			}
			LL[i][j] = (A[i][j] - sum)/LL[j][j];
		}
	}
	sum = 0.0;
	for(k=0;k<ReducedJDim-1;k++)
	{
		sum += LL[ReducedJDim-1][k]*LL[ReducedJDim-1][k];
	}
	LL[ReducedJDim-1][ReducedJDim-1] = sqrt(A[ReducedJDim-1][ReducedJDim-1] - sum);

	for(i=0;i<ReducedJDim-1;i++)
	{
		for(j=i+1;j<ReducedJDim;j++)
		{
			LL[i][j] = LL[j][i];
		}
	}
}

void Lyb2(double **L,double *y,double *b,int ReducedJDim)
{
	int i,j;
	y[0] = b[0]/L[0][0];
	for(i=1;i<ReducedJDim;i++)
	{
		y[i] = b[i];
		for(j=0;j<i;j++)
		{
			y[i] -= L[i][j]*y[j];
		}
		y[i] /= L[i][i];
	}
}

void ComputeNodalVariable(double **GlobMatrix,double *Nodal_Force,double *Nodal_Displacement,double *ReducedDisp,int *ReducedJRow,int JDim,int ReducedJDim)
{
	int i,j,k;
	double temp;
	for(i=0;i<ReducedJDim;i++)
	{
		j = ReducedJRow[i];
		Nodal_Displacement[j]=ReducedDisp[i];
	}
	for(i=0;i<JDim;i++)
	{
		temp = Nodal_Force[i];
		Nodal_Force[i] = 0.0;
		for(j=0;j<JDim;j++)
		{
			Nodal_Force[i] += GlobMatrix[i][j] * Nodal_Displacement[j];
		}
		Nodal_Force[i] -= temp;
	}
	for(i=0;i<ReducedJDim;i++)
	{
		j = ReducedJRow[i];
		Nodal_Force[j] = 0.0;
		for(k=0;k<JDim;k++)
		{
			Nodal_Force[j] += GlobMatrix[j][k] * Nodal_Displacement[k];
		}
	}
}
	
/* PostProcess */
void PostProcess(struct elem_structure *Elem,double ****ElemBIntPoint,double **JacobianIntPoint,double *Nodal_Displacement,double **D,double **StrainX,double **StrainY,double **StrainXY,double **StressX,double **StressY,double **StressXY,double **StressVonMises,double **StrainEnergy,double t,int ElemNum,int DoF,int NumIntPoint,int NumNodePerElem)
{
	int i,j,k;
	for(i=0;i<ElemNum;i++)
	{
		for(j=0;j<NumIntPoint*NumIntPoint;j++)
		{
			/* Initialize array */
			StrainX[i][j]  = 0.0;
			StrainY[i][j]  = 0.0;
			StrainXY[i][j] = 0.0;
			StressX[i][j]  = 0.0;
			StressY[i][j]  = 0.0;
			StressXY[i][j] = 0.0;
			StressVonMises[i][j] = 0.0;
			StrainEnergy[i][j]   = 0.0;
			/* Compute Strain = B * u Matrix */
			for(k=0;k<NumNodePerElem;k++)
			{
				StrainX[i][j]  += ElemBIntPoint[i][j][0][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].Node[k]] + ElemBIntPoint[i][j][0][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].Node[k]+1];
				StrainY[i][j]  += ElemBIntPoint[i][j][1][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].Node[k]] + ElemBIntPoint[i][j][1][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].Node[k]+1];
				StrainXY[i][j] += ElemBIntPoint[i][j][2][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].Node[k]] + ElemBIntPoint[i][j][2][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].Node[k]+1];
			}
			/* Compute Stress = D * Strain */
			StressX[i][j]  = D[0][0] * StrainX[i][j] + D[0][1] * StrainY[i][j] + D[0][2] * StrainXY[i][j];
			StressY[i][j]  = D[1][0] * StrainX[i][j] + D[1][1] * StrainY[i][j] + D[1][2] * StrainXY[i][j];
			StressXY[i][j] = D[2][0] * StrainX[i][j] + D[2][1] * StrainY[i][j] + D[2][2] * StrainXY[i][j];
			StressVonMises[i][j] = sqrt(pow(StressX[i][j],2)-StressX[i][j]*StressY[i][j]+pow(StressY[i][j],2)+3*pow(StressXY[i][j],2));
			StrainEnergy[i][j] = (StressX[i][j] * StrainX[i][j] + StressY[i][j] * StrainY[i][j] + StressXY[i][j] * StrainXY[i][j]) * t * JacobianIntPoint[i][j] / 2.0;
		}
	}
}
	
void PrintElemMatrix(FILE *fp3,double ***ElemMatrix,int KDim,int i)
{
	int j,k;
	fprintf(fp3,"Stiffness of Elem[%d]\n\n",i);
	for(j=0;j<KDim;j++)
	{
		for(k=0;k<KDim;k++)
		{
			fprintf(fp3," %20.12lf",ElemMatrix[i][j][k]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");
}

void PrintABD(FILE *fp3,double **AA,double **BB,double **DD,int DDim)
{
	int i,j;
	fprintf(fp3,"Matrix AA\n\n");
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			fprintf(fp3," %15.6lf",AA[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");

	fprintf(fp3,"Matrix BB\n\n");
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			fprintf(fp3," %15.6lf",BB[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");

	fprintf(fp3,"Matrix DD\n\n");
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			fprintf(fp3," %15.6lf",DD[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");
}

void PrintAC(FILE *fp3,double **A,double **C,int ADim,int DDim)
{
	int i,j;
	fprintf(fp3,"Matrix A\n\n");
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<ADim;j++)
		{
			fprintf(fp3," %15.6lf",A[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");

	fprintf(fp3,"Matrix C\n\n");
	for(i=0;i<DDim;i++)
	{
		for(j=0;j<DDim;j++)
		{
			fprintf(fp3," %15.6lf",C[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");
}

void PrintKe(FILE *fp3,double **Ke,int DimRow,int DimCol)
{
	int i,j;
	fprintf(fp3,"Matrix Ke\n\n");
	for(i=0;i<DimRow;i++)
	{
		for(j=0;j<DimCol;j++)
		{
			fprintf(fp3," %15.6lf",Ke[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fprintf(fp3,"\n");
}

void PrintArrary(FILE *fp3,double *Arrary,int Dim)
{
	int i,j;
	fprintf(fp3,"Arrary\n\n");
	for(i=0;i<Dim;i++)
	{
		fprintf(fp3," %15.6lf\n",Arrary[i]);
	}
	fprintf(fp3,"\n");
}

int main(int argc,char *argv[])
{
	int i,j,k;
	double temp1,temp2;
	FILE *fp1,*fp2,*fp3,*fp4;  /* fp1=.inp; fp2=.result; fp3=.sta; fp4=exsolu4node.dat */
	fp3 = fopen(argv[3],"w+"); 
	printf("Welcome to AME561_yilehu v1.0.\n\n");
	printf("You are now using the AME561_PS4node Solver.\n");
	printf("\n");
	/*********************************************************/
	/* Read Input File */
	/*********************************************************/
	int NumL,NumW;
	double L,W,r,a,t,E,v,Ta,Tb,ratioa;
	fp1 = fopen(argv[1],"r");
	fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf",&L,&W,&r,&a,&t,&E,&v,&Ta,&Tb,&NumL,&NumW,&ratioa);
	fclose(fp1);
	/*********************************************************/
	/* Discretization and connectivity */
	/*********************************************************/
	int ElemType = 2;
	int NodeNum,ElemNum;
	double dx,dy;
	struct node_structure *Node;
	struct elem_structure *Elem;

	NodeNum = (NumL+1)*(NumW+1);
	ElemNum = NumL*NumW;
	Node = (node_structure*)malloc(NodeNum*sizeof(node_structure));
	Elem = (elem_structure*)malloc(ElemNum*sizeof(elem_structure));
	Discretization(Node,L,W,NumL,NumW,dx,dy);
	Connectivity(Elem,NumL,NumW);
	/*********************************************************/
	/* Material Property */
	/*********************************************************/
	int DoF0,DoF,DDim,DDim2;
	double **D,**D2;
	double **AA,**BB,**DD;
	switch(ElemType)
	{
	case 1:
		DoF = 2;
		DDim = 3*DoF-3;
		D = (double**)malloc(DDim*sizeof(double*));
		for(j=0;j<DDim;j++)
			{*(D+j) = (double*)malloc(DDim*sizeof(double));}
		ElasticMatrixIsotropic(D,E,v);
		break;
	case 2:
		DoF0 = 2;
		DoF = 5;
		DDim = 3;
		DDim2 = 2;
		D = (double**)malloc(DDim*sizeof(double*));
		D2 = (double**)malloc(DDim2*sizeof(double*));
		AA = (double**)malloc(DDim*sizeof(double*));
		BB = (double**)malloc(DDim*sizeof(double*));
		DD = (double**)malloc(DDim*sizeof(double*));
		for(j=0;j<DDim;j++)
			{*(D+j) = (double*)malloc(DDim*sizeof(double));
			*(AA+j) = (double*)malloc(DDim*sizeof(double));
			*(BB+j) = (double*)malloc(DDim*sizeof(double));
			*(DD+j) = (double*)malloc(DDim*sizeof(double));}
		for(j=0;j<DDim2;j++)
			{*(D2+j) = (double*)malloc(DDim2*sizeof(double));}
		ElasticMatrixIsotropic(D,E,v);
		ComputeABDIsotropic(D,AA,BB,DD,DDim,t);
		PrintABD(fp3,AA,BB,DD,DDim);
		break;
	}
	/*********************************************************/
	/* Compute Elementary Stiffness Matrix */
	/*********************************************************/
	int NumIntPoint = 2;
	int NumNodePerElem = 4;
	int	KDim = DoF*NumNodePerElem;
	double **JacobianIntPoint1,**IntPointCoordX1,**IntPointCoordY1,**BTD1;
	double ***ElemMatrix1;
	double ****ElemBIntPoint1;

	double **A,**C;
	double **ElemMatrix2;

	switch(ElemType)
	{
	case 1:
		JacobianIntPoint1 = (double**)malloc(ElemNum*sizeof(double*));
		IntPointCoordX1 = (double**)malloc(ElemNum*sizeof(double*));
		IntPointCoordY1 = (double**)malloc(ElemNum*sizeof(double*));
		ElemMatrix1 = (double***)malloc(ElemNum*sizeof(double**));
		ElemBIntPoint1 = (double****)malloc(ElemNum*sizeof(double***));
		for(i=0;i<ElemNum;i++)
		{
			*(JacobianIntPoint1+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
			*(IntPointCoordX1+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
			*(IntPointCoordY1+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
			*(ElemMatrix1+i) = (double**)malloc(KDim*sizeof(double*));
			*(ElemBIntPoint1+i) = (double***)malloc(NumIntPoint*NumIntPoint*sizeof(double**));
			for(j=0;j<KDim;j++)
				{*(*(ElemMatrix1+i)+j) = (double*)malloc(KDim*sizeof(double));}
			for(j=0;j<NumIntPoint*NumIntPoint;j++)
				{*(*(ElemBIntPoint1+i)+j) = (double**)malloc(DDim*sizeof(double*));
				for(k=0;k<DDim;k++)
					{*(*(*(ElemBIntPoint1+i)+j)+k) = (double*)malloc(KDim*sizeof(double));}}
		}
		BTD1 = (double**)malloc(KDim*sizeof(double*));
		for(i=0;i<KDim;i++)
			{*(BTD1+i) = (double*)malloc(DDim*sizeof(double));}
		/* Compute ElemMatrix */
		for(i=0;i<ElemNum;i++)
		{
			InitializeElemMatrix(ElemBIntPoint1,JacobianIntPoint1,ElemMatrix1,DDim,KDim,NumIntPoint,i);
			BNumIntegral(Node,Elem,i,ElemBIntPoint1,JacobianIntPoint1,IntPointCoordX1,IntPointCoordY1);
			/* Assemble Elementary Stiffness Matrix*/
			for(j=0;j<NumIntPoint*NumIntPoint;j++)
			{
				InitializeBTDMatrix(BTD1,DDim,KDim);
				ComputeBTD(BTD1,ElemBIntPoint1,D,DDim,KDim,i,j);
				ComputeElemMatrix1(ElemMatrix1,BTD1,ElemBIntPoint1,JacobianIntPoint1,t,DDim,KDim,i,j);
			}
			PrintElemMatrix(fp3,ElemMatrix1,KDim,i);
		}
		break;
	case 2:
		A = (double**)malloc(DDim*sizeof(double*));
		C = (double**)malloc(DDim*sizeof(double*));
		ElemMatrix2 = (double**)malloc(KDim*sizeof(double*));
		for(i=0;i<DDim;i++)
		{
			*(A+i) = (double*)malloc(DoF0*DoF0*sizeof(double));
			*(C+i) = (double*)malloc(DDim*sizeof(double));
		}
		for(i=0;i<KDim;i++)
		{
			*(ElemMatrix2+i) = (double*)malloc(KDim*sizeof(double));
			for(j=0;j<KDim;j++){ElemMatrix2[i][j] = 0.0;}
		}
		ComputeAC(A,C,dx,dy);
		PrintAC(fp3,A,C,DoF0*DoF0,DDim);
		ComputeElemMatrix2(fp3,ElemMatrix2,AA,BB,DD,A,C,DDim,NumNodePerElem,DoF0,DoF,dx,dy);
		PrintKe(fp3,ElemMatrix2,KDim,KDim);
		break;
	}
	/*********************************************************/
	/* Assembly Global Stiffness Matrix */
	/*********************************************************/
	int JDim=DoF*NodeNum;
	double **GlobMatrix;
	GlobMatrix = (double**)malloc(JDim*sizeof(double*));
	for(i=0;i<JDim;i++)
	{
		*(GlobMatrix+i) = (double*)malloc(JDim*sizeof(double));
		for(j=0;j<JDim;j++){GlobMatrix[i][j] = 0.0;}
	}
	switch(ElemType)
	{
	case 1:
		ComputeGlobalMatrix1(GlobMatrix,ElemMatrix1,Elem,ElemNum,DoF,NumNodePerElem);
		break;
	case 2:
		ComputeGlobalMatrix2(GlobMatrix,ElemMatrix2,Elem,ElemNum,NodeNum,DoF0,DoF,NumNodePerElem);
		break;
	}
	PrintKe(fp3,GlobMatrix,JDim,JDim);
	/*********************************************************/
	/* Boundary Conditions */
	/*********************************************************/
	int *BCImplicitStatus;
	BCImplicitStatus = (int*)malloc(JDim*sizeof(int));
	switch(ElemType)
	{
	case 1:
		ComputeBC1(BCImplicitStatus,JDim,DoF,NumL,NumW);
		break;
	case 2:
		ComputeBC2(BCImplicitStatus,JDim,DoF0,DoF,NodeNum,NumL,NumW);
		break;
	}
	/*********************************************************/
	/* Compute Equivalent Nodal Force */
	/*********************************************************/
	double *Nodal_Force,*Nodal_Displacement;
	Nodal_Force = (double*)malloc(JDim*sizeof(double));
	Nodal_Displacement = (double*)malloc(JDim*sizeof(double));
	InitializeNodalVarible(Nodal_Force,Nodal_Displacement,DoF,NodeNum);
	switch(ElemType)
	{
	case 1:
		ComputeFeq1(Node,Nodal_Force,DoF,NumL,NumW,W,t,Ta,Tb);
		break;
	case 2:
		ComputeFeq2(Elem,Nodal_Force,DoF0,DoF,NodeNum,NumL,NumW,dx,dy,1.0);
		//i = (NodeNum-1)/2;
		//Nodal_Force[DoF0*NodeNum + (DoF-DoF0)*i+0] = 1000.0;
		break;
	}
	fp2 = fopen(argv[2],"w+");
	fprintf(fp2," Node               Coordx               Coordy                   Ux                   Uy                   Uz                   Rx                   Ry           Reaction x           Reaction y           Reaction z             Moment x             Moment y\n\n");
	for(i=0;i<NodeNum;i++)
	{
		fprintf(fp2,"%5d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,Node[i].X,Node[i].Y,Nodal_Displacement[DoF0*i+0],Nodal_Displacement[DoF0*i+1],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+0],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+1],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+2],Nodal_Force[DoF0*i+0],Nodal_Force[DoF0*i+1],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+0],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+1],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+2]);
	}
	fprintf(fp2,"\n");
	fclose(fp2);
	system("pause");
	/*********************************************************/
	/* Reduced Stiffness Matrix */
	/*********************************************************/
	int ReducedJDim;
	int *ReducedJRow;
	double **ReducedMatrix,**ReducedInvMatrix,*ReducedForce,*ReducedDisp,*ReducedArray;
	ReducedJDim = ComputeReducedJDim(BCImplicitStatus,JDim);
	ReducedJRow = (int*)malloc(ReducedJDim*sizeof(int));
	ReducedMatrix = (double**)malloc(ReducedJDim*sizeof(double*));
	ReducedInvMatrix = (double**)malloc(ReducedJDim*sizeof(double*));
	for(i=0;i<ReducedJDim;i++)
	{
		*(ReducedMatrix+i) = (double*)malloc(ReducedJDim*sizeof(double));
		*(ReducedInvMatrix+i)=(double*)malloc(ReducedJDim*sizeof(double));
	}
	ReducedForce = (double*)malloc(ReducedJDim*sizeof(double));
	ReducedDisp = (double*)malloc(ReducedJDim*sizeof(double));
	ReducedArray = (double*)malloc(ReducedJDim*sizeof(double));

	ComputeReducedJRow(BCImplicitStatus,ReducedJRow,JDim);
	ComputeReducedMatrix(GlobMatrix,Nodal_Force,ReducedMatrix,ReducedForce,ReducedDisp,ReducedJRow,ReducedJDim);
	PrintKe(fp3,ReducedMatrix,ReducedJDim,ReducedJDim);
	/*********************************************************/
	/* Solution */
	/*********************************************************/
	int START_CLOCK,END_CLOCK;
	double Running_Time;
	START_CLOCK = clock();
	//Inverse(ReducedMatrix,ReducedInvMatrix,ReducedForce,ReducedDisp,ReducedJDim);

	//GaussElimination(ReducedMatrix,ReducedInvMatrix,ReducedForce,ReducedJDim);
	//Uxy(ReducedInvMatrix,ReducedDisp,ReducedForce,ReducedJDim);

	//LUDecomposition(ReducedMatrix,ReducedInvMatrix,ReducedJDim);
	//Lyb(ReducedInvMatrix,ReducedArray,ReducedForce,ReducedJDim);
	//Uxy(ReducedInvMatrix,ReducedDisp,ReducedArray,ReducedJDim);

	LLDecomposition(ReducedMatrix,ReducedInvMatrix,ReducedJDim);
	Lyb2(ReducedInvMatrix,ReducedArray,ReducedForce,ReducedJDim);
	Uxy(ReducedInvMatrix,ReducedDisp,ReducedArray,ReducedJDim);

	END_CLOCK = clock();
	Running_Time = (double)(END_CLOCK - START_CLOCK)/CLOCKS_PER_SEC;
	
	ComputeNodalVariable(GlobMatrix,Nodal_Force,Nodal_Displacement,ReducedDisp,ReducedJRow,JDim,ReducedJDim);
	PrintKe(fp3,ReducedInvMatrix,ReducedJDim,ReducedJDim);
	PrintArrary(fp3,ReducedForce,ReducedJDim);
	PrintArrary(fp3,ReducedDisp,ReducedJDim);
	/*********************************************************/
	/* Postprocess */
	/*********************************************************/
	double **StrainX,**StrainY,**StrainXY;
	double **StressX,**StressY,**StressXY,**StressVonMises,**StrainEnergy;

	StrainX  = (double**)malloc(ElemNum*sizeof(double*));
	StrainY  = (double**)malloc(ElemNum*sizeof(double*));
	StrainXY = (double**)malloc(ElemNum*sizeof(double*));
	StressX  = (double**)malloc(ElemNum*sizeof(double*));
	StressY  = (double**)malloc(ElemNum*sizeof(double*));
	StressXY = (double**)malloc(ElemNum*sizeof(double*));
	StressVonMises = (double**)malloc(ElemNum*sizeof(double*));
	StrainEnergy   = (double**)malloc(ElemNum*sizeof(double*));
	for(i=0;i<ElemNum;i++)
	{
		*(StrainX+i)  = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StrainY+i)  = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StrainXY+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StressX+i)  = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StressY+i)  = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StressXY+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StressVonMises+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(StrainEnergy+i)   = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
	}
	//PostProcess(Elem,ElemBIntPoint1,JacobianIntPoint1,Nodal_Displacement,D,StrainX,StrainY,StrainXY,StressX,StressY,StressXY,StressVonMises,StrainEnergy,t,ElemNum,DoF,NumIntPoint,NumNodePerElem);
	/*********************************************************/
	printf("Geometry parameters:\n\n");
	printf("    L = %.4lf mm, W = %.4lf mm, r = %.4lf mm, a = %.4lf mm, t = %.4lf mm\n\n",L,W,r,a,t);
	printf("Material Property:\n\n");
	printf("    Young's Modulus E = %.4lf MPa, Poisson's Ratio v = %.4lf\n\n",E,v);
	printf("Boundary Condition:\n\n");
	printf("    Ta = %.4lf MPa, Tb = %.4lf MPa\n\n",Ta,Tb);
	printf("Mesh Density:\n\n");
	printf("    numr = %d, numa = %d, ratioa = %lf\n\n",NumL,NumW,ratioa);
	printf("Done!\n\n");
	printf("Node number  Element number  JDim            ReducedJDim     Max StressX     Total strain energy  Running_Time\n");
	printf("%-12d %-15d %-15d %-15d %-15e %-20e %-15e\n",NodeNum,ElemNum,JDim,ReducedJDim,temp1,temp2,Running_Time);
	/*********************************************************/
	/* Output Results */
	/*********************************************************/
	switch(ElemType)
	{
	case 1:
		temp1 = 0.0;
		temp2 = 0.0;
		fp2 = fopen(argv[2],"w+");
		fprintf(fp2," Node               Coordx               Coordy                   Ux                   Uy           Reaction x           Reaction y\n\n");
		for(i=0;i<NodeNum;i++)
		{
			fprintf(fp2,"%5d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,Node[i].X,Node[i].Y,Nodal_Displacement[DoF*i+0],Nodal_Displacement[DoF*i+1],Nodal_Force[DoF*i+0],Nodal_Force[DoF*i+1]);
		}
		fprintf(fp2,"\n");
		fprintf(fp2," Elem   IntPoint               Coordx               Coordy                  E11                  E22                  E12            SVonMises                  S11                  S22                  S12                  SEN\n\n");
		for(i=0;i<ElemNum;i++)
		{
			fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,1,IntPointCoordX1[i][0],IntPointCoordY1[i][0],StrainX[i][0],StrainY[i][0],StrainXY[i][0],StressVonMises[i][0],StressX[i][0],StressY[i][0],StressXY[i][0],StrainEnergy[i][0]);
			fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,2,IntPointCoordX1[i][2],IntPointCoordY1[i][2],StrainX[i][2],StrainY[i][2],StrainXY[i][2],StressVonMises[i][2],StressX[i][2],StressY[i][2],StressXY[i][2],StrainEnergy[i][2]);
			fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,3,IntPointCoordX1[i][1],IntPointCoordY1[i][1],StrainX[i][1],StrainY[i][1],StrainXY[i][1],StressVonMises[i][1],StressX[i][1],StressY[i][1],StressXY[i][1],StrainEnergy[i][1]);
			fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,4,IntPointCoordX1[i][3],IntPointCoordY1[i][3],StrainX[i][3],StrainY[i][3],StrainXY[i][3],StressVonMises[i][3],StressX[i][3],StressY[i][3],StressXY[i][3],StrainEnergy[i][3]);
			for(j=0;j<4;j++)
			{
				if (StressX[i][j] > temp1)
				{
					temp1 = StressX[i][j];
				}
				temp2 += StrainEnergy[i][j];
			}
		}
		fclose(fp2);
		fclose(fp3);
		break;
	case 2:
		fp2 = fopen(argv[2],"w+");
		fprintf(fp2," Node               Coordx               Coordy                   Ux                   Uy                   Uz                   Rx                   Ry           Reaction x           Reaction y           Reaction z             Moment x             Moment y\n\n");
		for(i=0;i<NodeNum;i++)
		{
			fprintf(fp2,"%5d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,Node[i].X,Node[i].Y,Nodal_Displacement[DoF0*i+0],Nodal_Displacement[DoF0*i+1],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+0],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+1],Nodal_Displacement[DoF0*NodeNum+(DoF-DoF0)*i+2],Nodal_Force[DoF0*i+0],Nodal_Force[DoF0*i+1],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+0],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+1],Nodal_Force[DoF0*NodeNum+(DoF-DoF0)*i+2]);
		}
		fprintf(fp2,"\n");
		fclose(fp2);
		break;
	}
	fclose(fp3);
	system("pause");
	///*********************************************************/
	//fp4 = fopen(argv[4],"r");
	//int *boundaryelem;
	//double *coordY,*solution;
	//boundaryelem=(int*)malloc((numa+enum1)*sizeof(int));
	//coordY=(double*)malloc(2*(numa+enum1)*sizeof(double));
	//solution=(double*)malloc(2*(numa+enum1)*sizeof(double));
	//for(i=0;i<numa;i++)
	//{
	//	boundaryelem[i] = (numr-1)*numa + i;
	//}
	//for(j=0;j<enum1;i++,j++)
	//{
	//	boundaryelem[i] = numr*numa+numr*enum1/2+(numr/2-1)*enum1 + j;
	//}
	//for(j=0;j<numa+enum1;j++)
	//{
	//	i=boundaryelem[j];
	//	coordY[2*j]=IntPointCoordY[i][1];
	//	coordY[2*j+1]=IntPointCoordY[i][3];
	//	solution[2*j]=StressX[i][1];
	//	solution[2*j+1]=StressX[i][3];
	//	//fprintf(fp3,"%5d %20.12lf %20.12lf %20.12lf\n",i,IntPointCoordX[i][1],IntPointCoordY[i][1],StressX[i][1]);
	//	//fprintf(fp3,"%5d %20.12lf %20.12lf %20.12lf\n",i,IntPointCoordX[i][3],IntPointCoordY[i][3],StressX[i][3]);
	//}
	////
	//int *elemlabel;
	//double *exX,*exY,*exStressX,*apStressX;
	//elemlabel = (int*)malloc(100*sizeof(int));
	//exX = (double*)malloc(100*sizeof(double));
	//exY = (double*)malloc(100*sizeof(double));
	//exStressX = (double*)malloc(100*sizeof(double));
	//apStressX = (double*)malloc(100*sizeof(double));
	//for(i=0;i<100;i++)
	//{
	//	fscanf(fp4,"%d %lf %lf %lf\n",&elemlabel[i],&exX[i],&exY[i],&exStressX[i]);
	//}
	//fclose(fp4);
	//printf("%lf\n",exStressX[99]);
	//for(i=0;i<100;i++)
	//{
	//	for(j=0;j<2*(numa+enum1)-1;j++)
	//	{
	//		if(coordY[j]-exY[i]<=0 && coordY[j+1]-exY[i]>=0)
	//		{
	//			apStressX[i]=solution[j]+(solution[j+1]-solution[j])*(exY[i]-coordY[j])/(coordY[j+1]-coordY[j]);
	//		}
	//		else if(coordY[0]>=exY[i])
	//		{
	//			apStressX[i]=solution[0]+(solution[0]-solution[1])*(exY[i]-coordY[0])/(coordY[0]-coordY[1]);
	//		}
	//		else if(coordY[2*(numa+enum1)-1]<=exY[i])
	//		{
	//			apStressX[i]=solution[2*(numa+enum1)-1]+(solution[2*(numa+enum1)-1]-solution[2*(numa+enum1)-2])*(exY[i]-coordY[2*(numa+enum1)-1])/(coordY[2*(numa+enum1)-1]-coordY[2*(numa+enum1)-2]);
	//		}
	//	}
	//	fprintf(fp3,"%5d %20.12lf\n",i,apStressX[i]);
	//}
	//fclose(fp3);
	///*********************************************************/
	//system("pause");
	return 0;
}