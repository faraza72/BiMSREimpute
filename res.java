vimport java.io.*;
import java.util.*;
public class res
{
	static int m=0,n=0;
	double[][] matrix;
	double[][] readMatrix()
	{
		
		try
        {
			BufferedReader in = new BufferedReader(new FileReader("missing2.txt"));	//reading files in specified directory
			String line, line1;
			if((line = in.readLine()) != null)	//file reading
			{
				String[] values = line.split(" ");
				m=Integer.parseInt(values[0]);
				n=Integer.parseInt(values[1]);
			}
			matrix = new double[m][n];
			
			for(int i = 0; i <m; i++)
			{
				line = in.readLine();
				String[] values = line.split(" ");
				for (int j = 0; j < n; j++)
				{    
					matrix[i][j]=Double.parseDouble(values[j]);
				}  
			}
			in.close();
		}
		catch( IOException ioException ) {}
		return matrix;
	}
	void display()
	{
		for(int i=0; i<m;i++)
		{
			for(int j=0; j<n; j++)
			{
				System.out.print(matrix[i][j]+" ");
			}
			System.out.println();
		}	
    }		
}