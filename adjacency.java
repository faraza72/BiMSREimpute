import java.io.*;
import java.util.*;
public class adjacency
{
	static int I,J;
	static int k=0;
	static int count1=0;
	private static double[][] G;
	int[][] iden_mat;
	int index, count;
	int[] miss_sample;
	adjacency miss_ob[];
	
	adjacency(){}
	
	adjacency(int index,int count,int miss_sample[])
	{
		
	this.index=index;
	this.count=count;
	this.miss_sample= new int[this.count];
	for(int i=0; i < this.count; i++)
		{
			this.miss_sample[i]=miss_sample[i];
		}
    }
	
	adjacency[] readIndex()
	{
		int a=0,b=0,index=0;
		
		int count=0;
		adjacency temp;
		res obj = new res();
        I = obj.m;
        J = obj.n;
		int miss_sample1[]=new int[J];
		//System.out.println(I+" "+J);
		try
        {
			BufferedReader ad = new BufferedReader(new FileReader("mised_data_info2(20%).txt"));
			String line;
			
				iden_mat = new int[I][J];
			
			for(int i = 0; i <I; i++)
			{
				for (int j = 0; j < J; j++)
				{    
					iden_mat[i][j]=1;
				}  
			}
			while((line = ad.readLine()) != null)
			{
				String[] iden = line.split("\t");
				a=Integer.parseInt(iden[0]);
				b=Integer.parseInt(iden[1]);
				//System.out.println(a+" "+b);
				iden_mat[a][b]=0;
			}
			
			for(int i = 0; i < I; i++)
			{
				for (int j = 0; j < J; j++)
				{    
					if(iden_mat[i][j]==0)
					{
						count1+=1;
						break;
					}
				} 
			}
			//System.out.println(count1);
			miss_ob = new adjacency[count1];
			
			int col_mis=0;
			for(int i=0; i<I; i++)
			{
				
				k=0;index=0;
				/*for(int j=0; j<J; j++)
					miss_sample1[j]=MAX; */
				for(int j=0; j<J; j++)
				{
					if(iden_mat[i][j]==0)
					{
						index=i;
						miss_sample1[k]=j;
						k++;
						count+=1;
					}
				}
				
				if(count>0)
				{
					miss_ob[col_mis]=new adjacency(index,count,miss_sample1);
					col_mis++;
				//	System.out.println(col_mis + "  " + i);
				}
				count=0;
			}
			//System.out.println(col_mis--);
			
			for(int i=0; i<col_mis; i++)
			{
				for(int j=0; j<col_mis-i-1; j++)
				{
					if(miss_ob[j].count>miss_ob[j+1].count)
					{
						temp=miss_ob[j];
						miss_ob[j]=miss_ob[j+1];
						miss_ob[j+1]=temp;
					}
				}
			//System.out.print(i+" " +miss_ob[i].index+"\t"+miss_ob[i].count+"\t");
				//System.out.print("  Column Missing: ");
				for(int j=0;j<miss_ob[i].count;j++)
				{
					//if(miss_ob[i].miss_sample[j]<MAX)
					//{
						//System.out.print(miss_ob[i].miss_sample[j] + " ");
					//}
				}
			//System.out.println();
				
			}
			ad.close();
		}
		
		catch( IOException ioException ) {}
		return miss_ob;
	}
}