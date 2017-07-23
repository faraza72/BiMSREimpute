public class Residue4
{

    static int I, J, K,mn=0;
    double zij, ziJ, zIj, zIJ, Hij = 0, Rij, R, sum=0, dist;
	double msdAdd=0,add=0, mean=0, deviation1=0, add1=0, add2=0, final_deviation=0;
    private static double[][] G, Z, F, S, Final;
    private static double[] H, H_selected, E, E_selected;
	private static int[] H_index, E_index;
	private static int[] selected_gene, selected_gene1;
	private static int[][] Rd;
	private static adjacency[] index_ob;
	public int target_col;
	int[] H_in;
	int mh=0, gene_no=4, missing_no=8;
	double H_diff,deviation=0;
	double[] H_bi;
	double[] msd= new double[missing_no];
	double[] orgValue=new double[missing_no];

    public static void main(String args[]) 
	{
		res obj = new res();
		adjacency obj1= new adjacency();
		Residue4 r1 = new Residue4();
        G = obj.readMatrix();
        I = obj.m;
        J = obj.n;
		H = new double[I-1];
        Z = new double[I][J];
		F = new double[2][J];
		Rd = new int[I][J]; 
		S = new double[I][J];
		index_ob = obj1.readIndex();
		K = obj1.count1;
		for(int i=0; i<I; i++)
		{
			for(int j=0; j<J; j++)
				Rd[i][j]= obj1.iden_mat[i][j];
		}
		r1.multi(G, Rd);
		new Residue4().calculateResidue();
    }
	
	void multi(double[][] G,int[][] Rd)
	{
		for(int i=0; i<I; i++)
		{
			for(int j=0; j<J; j++)
			{
				if(G[i][j] * Rd[i][j]==-0.0)
					S[i][j]= 0.0;
				else
					S[i][j]= G[i][j] * Rd[i][j];
			}
			//System.out.println();
		}
		//System.out.println(K);
		/*for(int i=0;i<K;i++)
			System.out.println(index[i]);*/
		/*for(int i=0; i<I; i++)
		{
			for(int j=0; j<J; j++)
			{
				System.out.print(S[i][j]+"\t");
			}
			System.out.println();
		}*/
	}

    void calculateResidue() 
	{
		Residue4 r1 = new Residue4();
		E = new double[I-1];
		selected_gene1 = new int[I];
		E_index = new int[I-1];
		E_selected = new double[I];
		selected_gene = new int[I];
		H_index = new int[I-1];
		H_selected = new double[I];
        int s,z=0,target,mis,count=0;
		double avg_res, avg=0, avg1=0;
		int [] final_selected;
		
		for(mis=0;mis<K;mis++)
		{
			count=0;
			z=0;
			mn=0;
			target=index_ob[mis].index;
			
			for (int q = 0; q < J; q++)
				F[z][q] = S[target][q];
			z++;
			
			for (int p = 0; p<I; p++) 
			{
				
				Hij=0;
				R=0;
				if (p != target) 
				{
					for (int q = 0; q < J; q++) 
					{
						F[z][q] = S[p][q];	
						//System.out.println(F[z][q]);
					}	
					
					for (int i = 0; i < 2; i++) 
					{
						for (int j = 0; j < J; j++) 
						{
							zij = F[i][j];
							//System.out.println("zij : " + zij);

							for (int k = 0; k < J; k++) 
							{
								ziJ += F[i][k];
							}
							ziJ /= J;
							// System.out.println("ziJ : " + ziJ);

							for (int k = 0; k < 2; k++) 
							{
								zIj += F[k][j];
							}
							zIj /= 2;
							//                System.out.println("zIj : " + zIj);
							for (int k = 0; k < 2; k++) 
							{
								for (int l = 0; l < J; l++) 
								{
									zIJ += F[k][l];
								}
							}
							zIJ = (zIJ / (2 * J));
							//System.out.println("zIJ : " + zIJ);
							
							Rij = zij - ziJ - zIj + zIJ;
							//System.out.println("Value : " + Rij);
							
							Z[i][j] = Rij;
							Rij=zij=zIj=ziJ=zIJ=0;
						}
					}
					
					/*System.out.println("Residue Matrix");
					for (int i = 0; i < I; i++) 
					{
						for (int j = 0; j < J; j++) 
						{
							System.out.print("\t" + Z[i][j]);
						}
						System.out.println();
					}*/
					
					for (int i = 0; i < 2; i++) 
					{
						for (int j = 0; j < J; j++) 
						{
							R = Z[i][j] * Z[i][j];
							Hij += R;
						}
					}
					Hij = Hij / (2 * J);
					H[mn]=Hij;
					H_index[mn]=p;
			
					//System.out.println("Calculating Euclidean Distance");
					
					for(int d=0;d<J;d++)
					{
						sum = sum + Math.pow((F[0][d]-F[1][d]),2.0);
					}
					dist=(Math.sqrt(sum))/J;
					E[mn]=dist;
					E_index[mn]=p;
					mn++;
					
					dist=0;
					sum=0;
				}
			}
			
			r1.sort(H, H_index, I-1);
			/*for(int i=0; i<gene_no; i++)
				System.out.println("Residue: "+H[i]+"\t"+H_index[i]); */
			r1.sort(E , E_index, I-1);
			/*for(int i=0; i<gene_no; i++)
				System.out.println("Euclidean: "+E[i]+"\t"+E_index[i]); */

			final_selected= new int[gene_no];
			
			Final= new double[gene_no+1][J];
			
			for(int i=0; i<gene_no; i++)
			{
				for(int j=0; j<gene_no; j++)
				{
					if(H_index[i]==E_index[j])
					{
						final_selected[i]=H_index[i];
						count++;
						break;
					}
				}
			}
			
			//System.out.println(count); 
			
			/*System.out.println("Final selected genes :");
			for(int i=0; i<gene_no; i++)
				System.out.print(final_selected[i]+"\t"); 
			System.out.println();*/
			
			//taking the target gene to the first row
			for(int x=0; x<J; x++)
				Final[0][x]=S[target][x];
			
			//selected genes
			for(int w=1; w<count+1; w++)
			{
				for(int x=0; x<J; x++)
				{
					Final[w][x]= S[final_selected[w-1]][x];
				}
			}
			
			/*System.out.println("Final");
			for(int i=0; i<count+1; i++)
			{
				for(int j=0; j<J; j++)
					System.out.print("\t"+Final[i][j]);
				System.out.println();
			}*/
			
			mean+=r1.ColResidue(Final, mis, target);	
			
			/*if(mis==0)
			{
				System.exit(0);
			}*/
		}
		r1.Cal_deviation(mean);
	}
	
	double ColResidue(double[][] Final_gene_no,int mis,int target)
	{
		Residue4 r1 = new Residue4();
		double[][] F_col=new double[gene_no+1][J];
		double[][] F_discarded=new double[gene_no+1][J];
		double[][] Z=new double[gene_no+1][J];
		double[][] F_new;
		double[] Euc_disc=new double[gene_no];
		double[] Euc_disc1;
		int[] index_track;
		double Rij=0,zij=0,zIj=0,ziJ=0,zIJ=0;
		double Hij=0,R=0,missing_value=0,dist1=0,dist2=0;
		double H_full=0,rowAvg=0,colAvg=0;
		double total_sum=0,sum1=0,sum2=0,total_sum1=0, mean1=0;
		int sp=0, sp1=0,z=0, y1=0;
		
		for(int i=0;i<index_ob[mis].count;i++)
		{
			
			target_col=index_ob[mis].miss_sample[i];
		//	System.out.println("target column :"+target_col);
			
			for (int q = 0; q < gene_no+1; q++) 
			{
                F_col[q][z] = Final_gene_no[q][target_col];
            }
			
			for (int p = 0; p < J; p++) 
			{
                if (p != target_col) 
				{
					z++;
					for (int q = 0; q < gene_no+1; q++) 
					{
						F_col[q][z] = Final_gene_no[q][p];
					}
				}
			}
			z=0;
			
			/*System.out.println("Final_gene_no...");
			for (int m=0;m<gene_no;m++)
			{
				 for(int j=0;j<J;j++)
				{
					System.out.print(Final_gene_no[m][j]+"\t");
				}
			    System.out.println();
			}*/
			
			/*System.out.println("F_col....");
			for (int m=0;m<gene_no+1;m++)
			{
				 for(int j=0;j<J;j++)
				{
					System.out.print(F_col[m][j]+"\t");
				}
			    System.out.println();
			}*/
			
			for (int i1=0;i1<gene_no+1;i1++)
			{
				 for(int j=0;j<J;j++)
				{
					F_discarded[i1][j]=F_col[i1][j];
				}
			}

			/*System.out.println("F_discarded....");
			for (int i1=0;i1<gene_no+1;i1++)
			{
				 for(int j=0;j<J;j++)
				{
					System.out.print(F_discarded[i1][j]+"\t");
				}
			    System.out.println();
			}*/
		
			//System.out.println(index_ob[mis].count);
			
			H_full= r1.RowResidue(F_discarded, J);
			
			int n=J-1;
			index_track=new int[J];
			
			for(int c=0;c<J;c++)
				index_track[c]=c;
			
			/*for(int c=0;c<n;c++)
				System.out.print(index_track[c]+"\t");
			System.out.println();*/
			
			while( n>=J/2)
			{
				
				H_bi= new double[n];
				H_in= new int[n];
				
				for(int x1=1; x1<n; x1++)
				{
					F_new= new double[gene_no+1][n]; 
					int k1=0;
					Rij=zij=zIj=ziJ=zIJ=Hij=R=0;
					
					for (int p = 0; p < n; p++) 
					{
						if (x1 != p ) 
						{
							for (int q = 0; q < gene_no+1; q++) 
							{
								F_new[q][k1] = F_discarded[q][p];
							}
							k1++;
						}
					}
					
					/*System.out.println("F_new for n: "+n);
					for (int l = 0; l < gene_no+1; l++)
					{ 
						for (int k = 0; k < n; k++) 
						{
							System.out.print(F_new[l][k]+"\t");
						}
						System.out.println();
					}*/
					
					for (int l = 0; l < gene_no+1; l++) 
					{
						
						for (int j = 0; j < n; j++) 
						{
							zij = F_new[l][j];
							//System.out.println("zij : " + zij);
							
							for (int k = 0; k < n; k++) 
							{
								ziJ += F_new[l][k];
							}
							ziJ /= n;
							// System.out.println("ziJ : " + ziJ);
							
							for (int k = 0; k < gene_no+1; k++) 
							{
								zIj += F_new[k][j];
							}
							zIj /= gene_no+1;
							//System.out.println("zIj : " + zIj);
							
							for (int k = 0; k < gene_no+1; k++) 
							{
								for (int r = 0; r < n; r++) 
								{
									zIJ += F_new[k][r];
								}
							}
							zIJ = (zIJ / ((gene_no+1)* n));
							//System.out.println("zIJ : " + zIJ);
							
							Rij = zij - ziJ - zIj + zIJ;
							//System.out.println("Value : " + Rij);
							
							Z[l][j] = Rij;
							//System.out.println("Zij :"+Z[l][j]);
							
							Rij = zij = zIj = ziJ = zIJ = 0;
						}
					}
					
					for (int w = 0; w < gene_no+1; w++) 
					{
						for (int x = 0; x < n; x++) 
						{
							R = Z[w][x] * Z[w][x];
							Hij += R;
						}
					}
					Hij = Hij / ((gene_no+1) * n);
					//System.out.println(Hij);
					
					 H_diff= H_full-Hij;
					
					if(H_diff >=0)
					{
						H_bi[y1]=H_diff;
						H_in[y1]=x1;
						//System.out.println("diff is:"+H_bi[y1]+"\t"+"index is :"+H_in[y1]);
						y1++;
					}
				}
				if(y1==0)
					break;
				r1.sort(H_bi, H_in, y1);	
				int f=H_in[y1-1];

				//System.out.println("deleted index = "+index_track[f]);
				/*for(int ab=0; ab<y1; ab++)
				System.out.println("diff is:"+H_bi[ab]+"\t"+"index is :"+H_in[ab]);
				System.out.println();*/
				
				//System.out.println("f = "+f);
				//System.out.println("Deleted column is: "+index_track[f]);
				
				int di=0;
				
				for(int c=0; c<gene_no+1; c++)
				{
					int dis_col=0;
					for(int c1=0; c1<n+1; c1++)
					{
						if(c1!=f)
						{
							F_discarded[c][dis_col]=F_discarded[c][c1];
							dis_col++;
						}
					}
				}
				
				for(int c=0;c<=n;c++)
				{
					if(c!=f)
					index_track[di++]=index_track[c];
				}
				
			/*	System.out.println("after deleting from index track...");
				for(int c=0;c<y1+1;c++)
				{
					System.out.print(index_track[c]+"\t");
				}
				System.out.println();
				
				//Final discarded matrix
				System.out.println("Final discarded for n: "+n);
				for (int l = 0; l < gene_no+1; l++)
				{ 
					for (int k = 0; k < n; k++) 
					{
						System.out.print(F_discarded[l][k]+"\t");
					}
					System.out.println();
				}
			*/	
				H_full= r1.RowResidue(F_discarded, n);
				n--;
				y1=0;	
			}

			// Calculating Row Average
			for(int row=1;row<gene_no+1;row++)
			{
				sum1=0;
				for(int col=0;col<n;col++)
				{
					sum1 =sum1 + Math.pow((F_discarded[0][col]-F_discarded[row][col]),2.0);
				}
				dist1=Math.sqrt(sum1);
			//	System.out.println("Sum Row: "+dist1);
				total_sum=total_sum+dist1;
				Euc_disc[sp]=dist1;
				sp++;
			}
			sp=0;
			
			//System.out.println(total_sum);
			
			for(int c=0;c<gene_no;c++)
			{
				rowAvg=rowAvg+((Euc_disc[c]/total_sum)*F_discarded[c+1][0]);
			}
			//System.out.println("rowAvg = "+rowAvg);
			
			//Calculating Column average
			Euc_disc1= new double[n];
			for(int col=1;col<n;col++)
			{
				sum2=0;
				for(int row=0;row<gene_no+1;row++)
				{
					sum2 = sum2 + Math.pow((F_discarded[row][0]-F_discarded[row][col]),2.0);
					
					//System.out.println(F_discarded[row][0] + "  " +F_discarded[row][col] );
				}
				//System.out.println(sum2 );
				dist2=Math.sqrt(sum2);
				//System.out.println("Sum Col: "+dist2);
				total_sum1=total_sum1+dist2;
				Euc_disc1[sp1]=dist2;
				sp1++;
			}
			sp1=0;
			//	System.out.println("total sum1 : "+total_sum1);
			
			for(int c=0;c<n;c++)
			{
				colAvg=colAvg+((Euc_disc1[c]/total_sum1)*F_discarded[0][c+1]);
			}
			//System.out.println("colAvg = "+colAvg);
			
			missing_value=(rowAvg + colAvg)/2;
			//System.out.println("missing_value = "+missing_value);
			
			mean1+=G[target][target_col];
			System.out.println("mean = "+mean1);
			
			System.out.println("mh = "+mh);
			
			orgValue[mh]=G[target][target_col];
			System.out.println("orgValue = "+orgValue[mh]);
			
			msd[mh]=Math.pow((G[target][target_col]-missing_value),2.0);
			System.out.println("msd["+mh+"] :"+msd[mh]);
			
			G[target][target_col]=missing_value;
			
			Rd[target][target_col]=1;
			
			
			/*for(int c=0;c<J;c++)
			{
				System.out.print(G[24][c]+"\t");
			}
			System.out.println();
			for(int c=0;c<J;c++)
			{
				System.out.print(Rd[24][c]+"\t");
			}
			System.out.println(); */
			mh++;
			r1.multi(G,Rd);
        }
		return mean1;
	}
	
	double RowResidue(double [][] F_discarded, int n)
	{
		//Final discarded matrix
			/*System.out.println("Final discarded");
			for (int l = 0; l < gene_no+1; l++)
			{ 
				for (int k = 0; k < n; k++) 
				{
					System.out.print(F_discarded[l][k]+"\t");
				}
				System.out.println();
			}
			*/
			double H_full;
			Hij = 0;
			R = 0;
			//System.out.println("\n\n Calculating Residue..\n");
			
			for (int l = 0; l < gene_no+1; l++) 
			{
				//System.out.println("5");
				for (int j = 0; j < n; j++) 
				{
					zij = F_discarded[l][j];
					//System.out.println("zij : " + zij);
					for (int k = 0; k < n; k++) 
					{
						ziJ += F_discarded[l][k];
					}
					ziJ /= n;
					// System.out.println("ziJ : " + ziJ);
					for (int k = 0; k < gene_no+1; k++) 
					{
						zIj += F_discarded[k][j];
					}
					zIj /= gene_no+1;
					//System.out.println("zIj : " + zIj);
					for (int k = 0; k < gene_no+1; k++) 
					{
						for (int r = 0; r < n; r++) 
						{
							zIJ += F_discarded[k][r];
						}
					}
					zIJ = (zIJ / ((gene_no+1)* n));
					//System.out.println("zIJ : " + zIJ);
					Rij = zij - ziJ - zIj + zIJ;
					//System.out.println("Value : " + Rij);
					Z[l][j] = Rij;
					//System.out.println("Zij :"+Z[l][j]);
					Rij = zij = zIj = ziJ = zIJ = 0;
				}
			}	                        
			for (int w = 0; w < gene_no+1; w++) 
			{
				for (int x = 0; x < n; x++) 
				{
					R = Z[w][x] * Z[w][x];
					Hij += R;
				}
			}
			Hij = Hij / ((gene_no+1) * n);
			//H_full= Hij;
			//H[p] = Hij;
			//System.out.println("Mean Squared residue of "+target+" and "+p+" is :" + H[p]);
			//System.out.println(" Residue for number of column = "+n+"\t" + Hij);
			
			return Hij;
	}
	
	void Cal_deviation(double mean)
	{
		System.out.println("Mean : "+mean);
		
		//System.out.println("Mh: "+mh);
		
		mean=mean/mh-1;
		//System.out.println("Mean : "+mean);
		for(int f=0;f<mh;f++)
		{
			deviation+=Math.pow((orgValue[f]-mean),2.0);
		}
		System.out.println("Deviation : "+deviation);
		deviation/=mh;
		deviation=Math.sqrt(deviation);
		
		for(int f=0;f<mh;f++)
		{
			msdAdd+=msd[f];
		}
		msdAdd/=mh;
		msdAdd=Math.sqrt(msdAdd);
		final_deviation=msdAdd/deviation;
		System.out.println("Final Deviation : "+final_deviation);
	}
	
	void sort(double[] H, int[] H_index, int row)
	{
		double temp;
		int temp1;
		for(int i=0; i<row; i++)
		{
			for(int j=0; j<row-i-1; j++)
			{
				if(H[j]>H[j+1])
				{
					temp=H[j];
					H[j]=H[j+1];
					H[j+1]=temp;
					temp1=H_index[j];
					H_index[j]=H_index[j+1];
					H_index[j+1]=temp1;
				}
			}
		}
	}
}