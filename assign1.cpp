//Arnav Jain
//13MA20011
//Operatioanl Research Lab assignment 1
//Solving Linear Programming Problem

#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

float d[10]={0};
float mat[10][10], b[10], temp[10][10];
float ans[10][10], z[10];
int R, C;

void solution(int n){
    float  c;
    int i, j, k;
    for(int i = 0 ; i < 10 ; i++)
    	d[i] = 0;
    // for(i = n - 1 ; i > 0 ; i--)     
    // {
    //     if(temp[i-1][0] < temp[i][0])
    //         for(j = 0 ; j <= n ; j++)
    //         {
    //             c = temp[i][j];
    //             temp[i][j] = temp[i - 1][j];
    //             temp[i - 1][j] = c;
    //         }
    // }
    
    //********* changing to upper triangular matrix*************//
    //********* Forward elimination process**************//
    for(k = 0 ; k < n - 1 ; k++)
        for(i = k ; i < n - 1 ; i++)
        {
            c= (temp[i + 1][k]/temp[k][k]) ;
            
            for(j = 0 ; j <= n ; j++)
                temp[i + 1][j] -= c*temp[k][j];
        }
    //***************** Backward Substitution method****************//
    for(i=0;i<n;i++)
    {
        for(j=0;j<=n;j++)
            printf("%6.1f",temp[i][j]);
        printf("\n");
    }
    for(i = n - 1 ; i >= 0 ; i--)
    {
        c = 0;
        for( j = i ; j <= n - 1 ; j++)
            c = c + temp[i][j]*d[j];
        d[i]=(temp[i][n] - c)/temp[i][i];
    }
}

void swapIt(int row1, int row2, int col)
{
	for (int i = 0; i < col; i++)
	{
		float temp1 = temp[row1][i];
		temp[row1][i] = temp[row2][i];
		temp[row2][i] = temp1;
	}
}


/* function for finding rank of matrix */
int rankOfMatrix(int R, int C)
{
	int rank = C;
	for (int row = 0; row < rank; row++)
	{
		if (temp[row][row])
		{
		for (int col = 0; col < R; col++)
		{
			if (col != row)
			{
				double mult = (double)temp[col][row] /
									temp[row][row];
				for (int i = 0; i < rank; i++)
				temp[col][i] -= mult * temp[row][i];
			}
		}
		}
		else
		{
			bool reduce = true;
			for (int i = row + 1; i < R; i++)
			{
				if (temp[i][row])
				{
					swapIt(row, i, rank);
					reduce = false;
					break ;
				}
			}
			if (reduce)
			{
				rank--;
				for (int i = 0; i < R; i ++)
					temp[i][row] = temp[i][rank];
			}
			row--;
		}
	}
	return rank;
}

void getMinimum(int solInd, int var){
    float mnanswer = 999999.0, curr = 0;
    for(int k = 0 ; k < solInd ; k++){
        int check = 1;
        curr = 0;
        for(int l = 0 ; l < var ; l++){
            if(ans[k][l] < 0){
                check = 0;
                break;
            }
        }
        if(check != 1)
            continue;
        for(int l = 0 ; l < var ; l++){
            curr = curr + ans[k][l]*z[l];
        }
        if(mnanswer > curr)
            mnanswer = curr;
        //find the answer here and print it here for each equation.
    }
    printf("The minimum is %6.1f\n", mnanswer);
}

void getMaximum(int solInd, int var){
    float mxanswer = 0.0, curr = 0;
    for(int k = 0 ; k < solInd ; k++){
        int check = 1;
        curr = 0;
        for(int l = 0 ; l < var ; l++){
            if(ans[k][l] < 0){
                check = 0;
                break;
            }
        }
        if(check != 1)
            continue;
        for(int l = 0 ; l < var ; l++){
            curr = curr + ans[k][l]*z[l];
        }
        if(mxanswer < curr)
            mxanswer = curr;
        //find the answer here and print it here for each equation.
    }
    printf("The maximum is %6.1f\n", mxanswer);
}
int main()
{
    int rank1, rank2, i, j, var, eqn;
    printf("Enter no of variables\n");
    scanf("%d",&var);
    
    printf("ENter no. of equations\n");
    scanf("%d",&eqn);
    //printf("Enter coefficient matrix - A \n");
    for(i = 0 ; i < eqn ; i++)
    {
    	printf("Enter coefficients with b of equation no %d\n" , i + 1);
        for(j = 0 ; j <= var ; j++)
        {
            scanf("%f",&mat[i][j]);
        }
    }
    //print_matrix(eqn,var+1);
    for(i = 0 ; i < eqn ; i++)
    {
        for(j = 0 ; j <= var ; j++){
        	temp[i][j] = mat[i][j];
        }
    }
    rank1 = rankOfMatrix(eqn , var);  
    for(i = 0 ; i < eqn ; i++)
    {
        for(j = 0 ; j <= var ; j++){
        	temp[i][j] = mat[i][j];
        }
    }
    rank2 = rankOfMatrix(eqn,var + 1);
    printf("Rank of A matrix is : %d \n", rank1);
    printf("Rank of A/b matrix is : %d \n", rank2);
    if(rank1 != rank2){
    	printf("Cant solve\n");
    	return 0;
    }
    int m = var - rank1;
    int pw2 = (int)pow(2, var);
    int cnt1 = 0, v, solInd = 0, ind = 0;
    for(int i = 1 ; i < pw2 ; i++){
    	v = i;
    	cnt1 = 0;
    	while(v != 0){
    		cnt1 = cnt1 + v%2;
    		v = v/2;
    	}
    	if(cnt1 != (var - m))
    		continue;
    	//cout << i << endl;
    	for(int k = 0 ; k < rank1 ; k++)
    	{
    		v = i;
    		ind = 0;
        	for(j = 0 ; j < var ; j++){
        		if(v%2 == 1){
        			temp[k][ind++] = mat[k][j];
        			//cout << temp[k][ind - 1] << " ";
        		}
        		v = v/2;
        	}
        	temp[k][ind] = mat[k][var];
        	//cout << temp[k][ind] << endl;
        }
        if(rankOfMatrix(cnt1, cnt1) != cnt1)
            continue;
        for(int k = 0 ; k < rank1 ; k++)
        {
            v = i;
            ind = 0;
            for(j = 0 ; j < var ; j++){
                if(v%2 == 1){
                    temp[k][ind++] = mat[k][j];
                    //cout << temp[k][ind - 1] << " ";
                }
                v = v/2;
            }
            temp[k][ind] = mat[k][var];
            //cout << temp[k][ind] << endl;
        }
        ind = 0;
        solution(var - m);
        v = i;
        
    	for(j = 0 ; j < var ; j++){
    		if(v%2 == 1){
    			//cout << d[ind] << " ";
    			ans[solInd][j] = d[ind++];
    		}
    		else{
    			ans[solInd][j] = 0;
    		}
    		v = v/2;            
    	}
    	solInd++;
    }
    //print all BFS
    for(int k = 0 ; k < solInd ; k++){
        int check = 1;
        for(int l = 0 ; l < var ; l++){
            if(ans[k][l] < 0){
                check = 0;
                break;
            }
        }
        if(check != 1)
            continue;
        printf("A basic feasible solution is ");
        for(int l = 0 ; l < var ; l++){
            printf("%6.2f ", ans[k][l]);
        }
        printf("\n");
        
    }
    // all answers are in array ans[][]
    printf("Enter coefficients of each variable as in objective function\n");
    for(i = 0 ; i < var ; i++){
    	//cout << "here" << endl;
        scanf("%f",&z[i]);
	}
    printf("Enter 1 to get maximum, 2 to get minimum, 3 to get both \n");
    int bl = 0;
    scanf("%d", &bl);
    if(bl&(1 << 1))
        getMinimum(solInd, var);
	if(bl&1)
        getMaximum(solInd, var);
    return 0;
}


