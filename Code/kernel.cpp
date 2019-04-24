#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <string.h>

const double Pi = 3.141592653589793;

struct node_structure
{
	double coordx;
	double coordy;
};

struct elem_structure
{
	int node[4];
};

void Discretization(struct node_structure *Node,double L,double W,int numL,int numW,double &dx,double &dy)
{
	int i,j,k;
	dx = L/(double)numL;
	dy = W/(double)numW;
	for(j=0;j<=numW;j++)
	{
		for(i=0;i<=numL;i++)
		{
			k = i + j*(numL+1);
			Node[k].coordx = (double)i * dx;
			Node[k].coordy = (double)j * dy;
		}
	}
}

void Connectivity(struct elem_structure *Elem,int numL,int numW)
{
	int i,j,k,m,n;
	for(j=0;j<numW;j++)
	{
		for(i=0;i<numL;i++)
		{
			k = i +     j*numL;
			m = i +     j*(numL+1);
			n = i + (j+1)*(numL+1);
			Elem[k].node[0] = m;
			Elem[k].node[1] = m+1;
			Elem[k].node[2] = n+1;
			Elem[k].node[3] = n;
		}
	}
}

void ElasticMatrixIsotropic(double **D,int DoF,double E,double v)
{
	double temp1;
	if(DoF == 2)
	{
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
}

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

void InitializeNodalVarible(double *Nodal_Force,double *Nodal_Displacement,int DoF,int nodenum)
{
	int i;
	for(i=0;i<DoF*nodenum;i++)
	{
		Nodal_Force[i] = 0.0;
		Nodal_Displacement[i] = 0.0;
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
			J[i][2*j+k] = (((Node[Elem[i].node[2]].coordx-Node[Elem[i].node[0]].coordx)*(Node[Elem[i].node[3]].coordy-Node[Elem[i].node[1]].coordy)+(Node[Elem[i].node[1]].coordx-Node[Elem[i].node[3]].coordx)*(Node[Elem[i].node[2]].coordy-Node[Elem[i].node[0]].coordy))+((Node[Elem[i].node[2]].coordx-Node[Elem[i].node[3]].coordx)*(Node[Elem[i].node[0]].coordy-Node[Elem[i].node[1]].coordy)+(Node[Elem[i].node[1]].coordx-Node[Elem[i].node[0]].coordx)*(Node[Elem[i].node[2]].coordy-Node[Elem[i].node[3]].coordy))*kexi[j] + ((Node[Elem[i].node[2]].coordx-Node[Elem[i].node[1]].coordx)*(Node[Elem[i].node[3]].coordy-Node[Elem[i].node[0]].coordy)+(Node[Elem[i].node[0]].coordx-Node[Elem[i].node[3]].coordx)*(Node[Elem[i].node[2]].coordy-Node[Elem[i].node[1]].coordy))*eta[k]) / 8.0;
			X[i][2*j+k] = 0.25*(1.0-kexi[j])*(1.0-eta[k])*Node[Elem[i].node[0]].coordx + 0.25*(1.0+kexi[j])*(1.0-eta[k])*Node[Elem[i].node[1]].coordx + 0.25*(1.0+kexi[j])*(1.0+eta[k])*Node[Elem[i].node[2]].coordx + 0.25*(1.0-kexi[j])*(1.0+eta[k])*Node[Elem[i].node[3]].coordx;
			Y[i][2*j+k] = 0.25*(1.0-kexi[j])*(1.0-eta[k])*Node[Elem[i].node[0]].coordy + 0.25*(1.0+kexi[j])*(1.0-eta[k])*Node[Elem[i].node[1]].coordy + 0.25*(1.0+kexi[j])*(1.0+eta[k])*Node[Elem[i].node[2]].coordy + 0.25*(1.0-kexi[j])*(1.0+eta[k])*Node[Elem[i].node[3]].coordy;
			B[i][2*j+k][0][0] = ((Node[Elem[i].node[1]].coordy-Node[Elem[i].node[3]].coordy)+(Node[Elem[i].node[3]].coordy-Node[Elem[i].node[2]].coordy)*kexi[j]+(Node[Elem[i].node[2]].coordy-Node[Elem[i].node[1]].coordy)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][2] = ((Node[Elem[i].node[2]].coordy-Node[Elem[i].node[0]].coordy)+(Node[Elem[i].node[2]].coordy-Node[Elem[i].node[3]].coordy)*kexi[j]+(Node[Elem[i].node[0]].coordy-Node[Elem[i].node[3]].coordy)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][4] = ((Node[Elem[i].node[3]].coordy-Node[Elem[i].node[1]].coordy)+(Node[Elem[i].node[0]].coordy-Node[Elem[i].node[1]].coordy)*kexi[j]+(Node[Elem[i].node[3]].coordy-Node[Elem[i].node[0]].coordy)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][0][6] = ((Node[Elem[i].node[0]].coordy-Node[Elem[i].node[2]].coordy)+(Node[Elem[i].node[1]].coordy-Node[Elem[i].node[0]].coordy)*kexi[j]+(Node[Elem[i].node[1]].coordy-Node[Elem[i].node[2]].coordy)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][1] = ((Node[Elem[i].node[3]].coordx-Node[Elem[i].node[1]].coordx)+(Node[Elem[i].node[2]].coordx-Node[Elem[i].node[3]].coordx)*kexi[j]+(Node[Elem[i].node[1]].coordx-Node[Elem[i].node[2]].coordx)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][3] = ((Node[Elem[i].node[0]].coordx-Node[Elem[i].node[2]].coordx)+(Node[Elem[i].node[3]].coordx-Node[Elem[i].node[2]].coordx)*kexi[j]+(Node[Elem[i].node[3]].coordx-Node[Elem[i].node[0]].coordx)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][5] = ((Node[Elem[i].node[1]].coordx-Node[Elem[i].node[3]].coordx)+(Node[Elem[i].node[1]].coordx-Node[Elem[i].node[0]].coordx)*kexi[j]+(Node[Elem[i].node[0]].coordx-Node[Elem[i].node[3]].coordx)*eta[k]) / J[i][2*j+k] / 8.0;
			B[i][2*j+k][1][7] = ((Node[Elem[i].node[2]].coordx-Node[Elem[i].node[0]].coordx)+(Node[Elem[i].node[0]].coordx-Node[Elem[i].node[1]].coordx)*kexi[j]+(Node[Elem[i].node[2]].coordx-Node[Elem[i].node[1]].coordx)*eta[k]) / J[i][2*j+k] / 8.0;
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
			
void ComputeElemMatrix(double ***ElemMatrix,double **BTD,double ****ElemBIntPoint,double **JacobianIntPoint,double t,int DDim,int KDim,int i,int j)
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

void ComputeGlobalMatrix(double **GlobMatrix,double ***ElemMatrix,struct elem_structure *Elem,int elemnum,int DoF,int NumNodePerElem)
{
	int i,j,k;
	for(i=0;i<elemnum;i++)
	{
		for(j=0;j<NumNodePerElem;j++)
		{
			for(k=0;k<NumNodePerElem;k++)
			{
				GlobMatrix[DoF*Elem[i].node[j]+0][DoF*Elem[i].node[k]+0] += ElemMatrix[i][DoF*j+0][DoF*k+0];
				GlobMatrix[DoF*Elem[i].node[j]+0][DoF*Elem[i].node[k]+1] += ElemMatrix[i][DoF*j+0][DoF*k+1];
				GlobMatrix[DoF*Elem[i].node[j]+1][DoF*Elem[i].node[k]+0] += ElemMatrix[i][DoF*j+1][DoF*k+0];
				GlobMatrix[DoF*Elem[i].node[j]+1][DoF*Elem[i].node[k]+1] += ElemMatrix[i][DoF*j+1][DoF*k+1];
			}
		}
	}
}

void ComputeFeq(struct node_structure *Node,double *Nodal_Force,int DoF,int numL,int numW,double W,double t,double Ta,double Tb)
{
	int i,j,k;
	double *T;
	T = (double*)malloc((numW+1)*sizeof(double));
	for(i=0;i<numW+1;i++)
	{
		j = (i+1)*(numL+1) - 1;
		T[i] = Tb + Node[j].coordy * (Ta - Tb) / W;
	}
	for(i=0;i<numW;i++)
	{
		j = (i+1)*(numL+1) - 1;
		k = (i+2)*(numL+1) - 1;
		//Nodal_Force[DoF*j+0] += t*(Node[k].coordy-Node[j].coordy) * (2*T[i+0] + T[i+1])/6.0;
		//Nodal_Force[DoF*k+0] += t*(Node[k].coordy-Node[j].coordy) * (T[i+0] + 2*T[i+1])/6.0;
		Nodal_Force[DoF*j+0] += t*(Node[k].coordy-Node[j].coordy) * (T[i+0] + T[i+1])/4.0;
		Nodal_Force[DoF*k+0] += t*(Node[k].coordy-Node[j].coordy) * (T[i+1] + T[i+0])/4.0;
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

void ComputeReducedJRow(int *BCImplicitStatus,int *ReducedJRow,int JDim,int DoF,int numL,int numW)
{
	int i,j;
	for(i=0;i<JDim;i++)
	{
		BCImplicitStatus[i] = 0;
	}
	for(i=0;i<numW+1;i++)
	{
		j = i*(numL+1);
		BCImplicitStatus[DoF*j+0] = 1;
		BCImplicitStatus[DoF*j+1] = 1;
	}
	for(i=0,j=0;i<JDim;i++)
	{
		if(BCImplicitStatus[i] == 0)
		{
			ReducedJRow[j] = i;
			j++;
		}
	}
}

void ComputeNodalVariable(double **GlobMatrix,double *Nodal_Force,double *Nodal_Displacement,double **ReducedInvMatrix,double *ReducedForce,double *ReducedDisp,int *ReducedJRow,int JDim,int ReducedJDim)
{
	int i,j;
	for(i=0;i<ReducedJDim;i++)
	{
		for(j=0;j<ReducedJDim;j++)
		{
			ReducedDisp[i] += ReducedInvMatrix[i][j] * ReducedForce[j];
		}
	}
	for(i=0;i<ReducedJDim;i++)
	{
		j = ReducedJRow[i];
		Nodal_Displacement[j]=ReducedDisp[i];
	}
	for(i=0;i<JDim;i++)
	{
		for(j=0;j<JDim;j++)
		{
			Nodal_Force[i] += GlobMatrix[i][j] * Nodal_Displacement[j];
		}
	}
}
	
void Inverse(double **R,double **Rinv, int num)
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
}

void PostProcess(struct elem_structure *Elem,double ****ElemBIntPoint,double **JacobianIntPoint,double *Nodal_Displacement,double **D,double **StrainX,double **StrainY,double **StrainXY,double **StressX,double **StressY,double **StressXY,double **StressVonMises,double **StrainEnergy,double t,int elemnum,int DoF,int NumIntPoint,int NumNodePerElem)
{
	int i,j,k;
	for(i=0;i<elemnum;i++)
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
				StrainX[i][j]  += ElemBIntPoint[i][j][0][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].node[k]] + ElemBIntPoint[i][j][0][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].node[k]+1];
				StrainY[i][j]  += ElemBIntPoint[i][j][1][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].node[k]] + ElemBIntPoint[i][j][1][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].node[k]+1];
				StrainXY[i][j] += ElemBIntPoint[i][j][2][DoF*k+0]*Nodal_Displacement[DoF*Elem[i].node[k]] + ElemBIntPoint[i][j][2][DoF*k+1]*Nodal_Displacement[DoF*Elem[i].node[k]+1];
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

int main(int argc,char *argv[])
{
	int i,j,k,l,m,n;
	double temp1,temp2,temp3,temp4;
	FILE *fp1,*fp2,*fp3,*fp4;  /* fp1=.inp; fp2=.result; fp3=.sta; fp4=exsolu4node.dat */
	fp3 = fopen(argv[3],"w+"); 
	printf("Welcome to AME561_yilehu v1.0.\n\n");
	printf("You are now using the AME561_PS4node Solver.\n");
	printf("\n");
	/*********************************************************/
	/* Read Input File */
	/*********************************************************/
	int numL,numW;
	double L,W,r,a,t,E,v,Ta,Tb,ratioa;
	fp1 = fopen(argv[1],"r");
	fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf",&L,&W,&r,&a,&t,&E,&v,&Ta,&Tb,&numL,&numW,&ratioa);
	fclose(fp1);
	/*********************************************************/
	/* Discretization and connectivity */
	/*********************************************************/
	int nodenum,elemnum;
	struct node_structure *Node;
	struct elem_structure *Elem;
	double dx,dy;
	nodenum = (numL+1)*(numW+1);
	elemnum = numL*numW;
	Node = (node_structure*)malloc(nodenum*sizeof(node_structure));
	Elem = (elem_structure*)malloc(elemnum*sizeof(elem_structure));
	Discretization(Node,L,W,numL,numW,dx,dy);
	Connectivity(Elem,numL,numW);
	/*********************************************************/
	/* Material Property */
	/*********************************************************/
	int DoF = 2;
	int DDim = 3*DoF-3;
	double **D;
	D = (double**)malloc(DDim*sizeof(double*));
	for(j=0;j<DDim;j++)
		{*(D+j) = (double*)malloc(DDim*sizeof(double));}
	ElasticMatrixIsotropic(D,DoF,E,v);
	/*********************************************************/
	/* Compute Elementary Stiffness Matrix */
	/*********************************************************/
	int NumIntPoint = 2;
	int NumNodePerElem = 4;
	int KDim = DoF*NumNodePerElem;
	double **JacobianIntPoint,**IntPointCoordX,**IntPointCoordY,**BTD;
	double ***ElemMatrix;
	double ****ElemBIntPoint;
	JacobianIntPoint = (double**)malloc(elemnum*sizeof(double*));
	IntPointCoordX = (double**)malloc(elemnum*sizeof(double*));
	IntPointCoordY = (double**)malloc(elemnum*sizeof(double*));
	ElemMatrix = (double***)malloc(elemnum*sizeof(double**));
	ElemBIntPoint = (double****)malloc(elemnum*sizeof(double***));
	for(i=0;i<elemnum;i++)
	{
		*(JacobianIntPoint+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(IntPointCoordX+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(IntPointCoordY+i) = (double*)malloc(NumIntPoint*NumIntPoint*sizeof(double));
		*(ElemMatrix+i) = (double**)malloc(KDim*sizeof(double*));
		*(ElemBIntPoint+i) = (double***)malloc(NumIntPoint*NumIntPoint*sizeof(double**));
		for(j=0;j<KDim;j++)
			{*(*(ElemMatrix+i)+j) = (double*)malloc(KDim*sizeof(double));}
		for(j=0;j<NumIntPoint*NumIntPoint;j++)
			{*(*(ElemBIntPoint+i)+j) = (double**)malloc(DDim*sizeof(double*));
			for(k=0;k<DDim;k++)
				{*(*(*(ElemBIntPoint+i)+j)+k) = (double*)malloc(KDim*sizeof(double));}}
	}
	BTD = (double**)malloc(KDim*sizeof(double*));
	for(i=0;i<KDim;i++)
		{*(BTD+i) = (double*)malloc(DDim*sizeof(double));}
	/* Compute ElemMatrix */
	for(i=0;i<elemnum;i++)
	{
		InitializeElemMatrix(ElemBIntPoint,JacobianIntPoint,ElemMatrix,DDim,KDim,NumIntPoint,i);
		BNumIntegral(Node,Elem,i,ElemBIntPoint,JacobianIntPoint,IntPointCoordX,IntPointCoordY);
		/* Assemble Elementary Stiffness Matrix*/
		for(j=0;j<NumIntPoint*NumIntPoint;j++)
		{
			InitializeBTDMatrix(BTD,DDim,KDim);
			ComputeBTD(BTD,ElemBIntPoint,D,DDim,KDim,i,j);
			ComputeElemMatrix(ElemMatrix,BTD,ElemBIntPoint,JacobianIntPoint,t,DDim,KDim,i,j);
		}
		PrintElemMatrix(fp3,ElemMatrix,KDim,i);
	}
	/*********************************************************/
	/* Assembly Global Stiffness Matrix */
	/*********************************************************/
	int JDim=DoF*nodenum;
	double **GlobMatrix;
	GlobMatrix = (double**)malloc(JDim*sizeof(double*));
	for(i=0;i<JDim;i++)
	{
		*(GlobMatrix+i) = (double*)malloc(JDim*sizeof(double));
		for(j=0;j<JDim;j++){GlobMatrix[i][j] = 0.0;}
	}
	ComputeGlobalMatrix(GlobMatrix,ElemMatrix,Elem,elemnum,DoF,NumNodePerElem);
	/*********************************************************/
	/* Compute Equivalent Nodal Force */
	/*********************************************************/
	int *BCImplicitStatus;
	double *Nodal_Force,*Nodal_Displacement;
	BCImplicitStatus = (int*)malloc(JDim*sizeof(int));
	Nodal_Force = (double*)malloc(JDim*sizeof(double));
	Nodal_Displacement = (double*)malloc(JDim*sizeof(double));
	
	InitializeNodalVarible(Nodal_Force,Nodal_Displacement,DoF,nodenum);
	ComputeFeq(Node,Nodal_Force,DoF,numL,numW,W,t,Ta,Tb);
	/*********************************************************/
	/* Reduced Stiffness Matrix */
	/*********************************************************/
	int ReducedJDim = JDim - DoF*(numW+1);
	int *ReducedJRow;
	double **ReducedMatrix,**ReducedInvMatrix,*ReducedForce,*ReducedDisp;
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

	ComputeReducedJRow(BCImplicitStatus,ReducedJRow,JDim,DoF,numL,numW);
	ComputeReducedMatrix(GlobMatrix,Nodal_Force,ReducedMatrix,ReducedForce,ReducedDisp,ReducedJRow,ReducedJDim);
	/*********************************************************/
	/* Solution */
	/*********************************************************/
	Inverse(ReducedMatrix,ReducedInvMatrix,ReducedJDim);
	ComputeNodalVariable(GlobMatrix,Nodal_Force,Nodal_Displacement,ReducedInvMatrix,ReducedForce,ReducedDisp,ReducedJRow,JDim,ReducedJDim);
	/*********************************************************/
	/* Postprocess */
	/*********************************************************/
	double **StrainX,**StrainY,**StrainXY;
	double **StressX,**StressY,**StressXY,**StressVonMises,**StrainEnergy;
	StrainX  = (double**)malloc(elemnum*sizeof(double*));
	StrainY  = (double**)malloc(elemnum*sizeof(double*));
	StrainXY = (double**)malloc(elemnum*sizeof(double*));
	StressX  = (double**)malloc(elemnum*sizeof(double*));
	StressY  = (double**)malloc(elemnum*sizeof(double*));
	StressXY = (double**)malloc(elemnum*sizeof(double*));
	StressVonMises = (double**)malloc(elemnum*sizeof(double*));
	StrainEnergy   = (double**)malloc(elemnum*sizeof(double*));
	for(i=0;i<elemnum;i++)
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
	PostProcess(Elem,ElemBIntPoint,JacobianIntPoint,Nodal_Displacement,D,StrainX,StrainY,StrainXY,StressX,StressY,StressXY,StressVonMises,StrainEnergy,t,elemnum,DoF,NumIntPoint,NumNodePerElem);
	/*********************************************************/
	/* Output Results */
	/*********************************************************/
	temp1 = 0.0;
	temp2 = 0.0;
	fp2 = fopen(argv[2],"w+");
	fprintf(fp2," Node               Coordx               Coordy                   Ux                   Uy           Reaction x           Reaction y\n\n");
	for(i=0;i<nodenum;i++)
	{
		fprintf(fp2,"%5d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,Node[i].coordx,Node[i].coordy,Nodal_Displacement[DoF*i+0],Nodal_Displacement[DoF*i+1],Nodal_Force[DoF*i+0],Nodal_Force[DoF*i+1]);
	}
	fprintf(fp2,"\n");
	fprintf(fp2," Elem   IntPoint               Coordx               Coordy                  E11                  E22                  E12            SVonMises                  S11                  S22                  S12                  SEN\n\n");
	for(i=0;i<elemnum;i++)
	{
		fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,1,IntPointCoordX[i][0],IntPointCoordY[i][0],StrainX[i][0],StrainY[i][0],StrainXY[i][0],StressVonMises[i][0],StressX[i][0],StressY[i][0],StressXY[i][0],StrainEnergy[i][0]);
		fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,2,IntPointCoordX[i][2],IntPointCoordY[i][2],StrainX[i][2],StrainY[i][2],StrainXY[i][2],StressVonMises[i][2],StressX[i][2],StressY[i][2],StressXY[i][2],StrainEnergy[i][2]);
		fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,3,IntPointCoordX[i][1],IntPointCoordY[i][1],StrainX[i][1],StrainY[i][1],StrainXY[i][1],StressVonMises[i][1],StressX[i][1],StressY[i][1],StressXY[i][1],StrainEnergy[i][1]);
		fprintf(fp2,"%5d %10d %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",i,4,IntPointCoordX[i][3],IntPointCoordY[i][3],StrainX[i][3],StrainY[i][3],StrainXY[i][3],StressVonMises[i][3],StressX[i][3],StressY[i][3],StressXY[i][3],StrainEnergy[i][3]);
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
	/*********************************************************/
	printf("Geometry parameters:\n\n");
	printf("    L = %.4lf mm, W = %.4lf mm, r = %.4lf mm, a = %.4lf mm, t = %.4lf mm\n\n",L,W,r,a,t);
	printf("Material Property:\n\n");
	printf("    Young's Modulus E = %.4lf MPa, Poisson's Ratio v = %.4lf\n\n",E,v);
	printf("Boundary Condition:\n\n");
	printf("    Ta = %.4lf MPa, Tb = %.4lf MPa\n\n",Ta,Tb);
	printf("Mesh Density:\n\n");
	printf("    numr = %d, numa = %d, ratioa = %lf\n\n",numL,numW,ratioa);
	printf("Done!\n\n");
	printf("Node number  Element number  Max StressX     Total strain energy\n");
	printf("%-12d %-15d %-15e %-15e\n",nodenum,elemnum,temp1,temp2);
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