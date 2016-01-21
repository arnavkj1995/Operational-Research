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

typedef long long ll;
typedef pair <int,int> pii;
typedef vector <int> vi;

#define rep(i, n) for(int i = 0; i < (n); ++i)
#define forn(i, a, b) for(int i = (a); i < (b); ++i)
#define ford(i, a, b) for(int i = (a); i >= (b); --i)
#define fore(i, a, b) forn(i, a, b + 1)

#define pb push_back
#define mp make_pair
#define ff first
#define ss second
#define all(c) c.begin(), c.end()
#define fill(a, v) memset(a, v, sizeof(a))
#define sz(a) ((int)a.size())

#define gl(x) cin >> x
#define gi(x) scanf("%d", &x)
#define pls(x) cout << x << " "
#define pln(x) cout << x << "\n"
#define pis(x) printf("%d ", x)
#define pin(x) printf("%d\n", x)
#define pnl printf("\n")
#define dbn cerr << "\n"
#define dbg(x) cerr << #x << " : " << x << " "
#define dbs(x) cerr << x << " "

#define foreach(c, it) for(__typeof(c.begin()) it = c.begin(); it != c.end(); ++it)

float d[10]={0};
float mat[10][10], b[10], temp[10][10], constants[10];
float ans[10][10], z[10];
int R, C;
int countVal=0,countVal1=0;

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

void solution(int n){
    float  c;
    int i, j, k;
    for(int i = 0 ; i < 10 ; i++)
        d[i] = 0;
    for(k = 0 ; k < n - 1 ; k++)
        for(i = k ; i < n - 1 ; i++)
        {
            c= (temp[i + 1][k]/temp[k][k]) ;
            
            for(j = 0 ; j <= n ; j++)
                temp[i + 1][j] -= c*temp[k][j];
        }
    for(i = n - 1 ; i >= 0 ; i--)
    {
        c = 0;
        for( j = i ; j <= n - 1 ; j++)
            c = c + temp[i][j]*d[j];
        d[i]=(temp[i][n] - c)/temp[i][i];
    }
}

void getMinimum(int solInd, int var){
    float mnanswer = 999999.0, curr = 0, cntmin = 0;
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
        if(mnanswer > curr){
            mnanswer = curr;
            cntmin = 1;
        }
        else if(mnanswer == curr){
            cntmin++;
        }
    }
    if(cntmin == 1)
        printf("The minimum is %6.1f\n", mnanswer);
    else
        printf("We have infinite solution\n");
}

void getMaximum(int solInd, int var){
    float mxanswer = 0.0, curr = 0, cntmax = 0;
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

        if(mxanswer < curr){
            mxanswer = curr;
            cntmax = 1;
        }
        else if(mxanswer == curr)
            cntmax++;
    }
    if(cntmax == 1)
        printf("The maximum is %6.1f\n", mxanswer);
    else
        printf("We have infinite solution\n");
}

void reduce(int rank, int var, int eqn){
    int i, j, v, cnt1, ind = 0, pw = (int)pow(2, eqn + 1);
    for(i = 1 ; i < pw ; i++){
        v = i;
        cnt1 = 0;
        while(v != 0){
            cnt1 = cnt1 + v%2;
            v = v/2;
        }
        if(cnt1 != rank)
            continue;
        ind = 0;
        v = i;
        for(int k = 0 ; k < eqn ; k++)
        {            
            if(v%2 != 1){
                v = v/2;
                continue;
            }
            for(j = 0 ; j <= var ; j++){             
                temp[ind][j] = mat[k][j];
            }
            v = v/2;
            ind++;
        }
        if(rankOfMatrix(rank, var + 1) == rank){
            break;
        }
    }
    ind = 0;
    v = i;
    for(int k = 0 ; k < eqn ; k++)
    {       
        if(v%2 != 1){
            v = v/2;
            continue;
        }
        for(j = 0 ; j <= var ; j++){                
            temp[ind][j] = mat[k][j];
        }
        v = v/2;
        ind++;
    }
    for(int k = 0 ; k < rank ; k++)
    {
        for(j = 0 ; j <= var ; j++){                
            mat[k][j] = temp[k][j];
        }
    }
}

int main()
{
    int rank1, rank2, i, j, k, in_var, var, eqn;
    string inequality;  
    printf("Enter no of variables\n");
    scanf("%d",&in_var);

    var= in_var; 
    
    printf("Enter no. of equations\n");
    scanf("%d",&eqn);

    for(i = 0 ; i < eqn ; i++)
    {
        printf("Enter coefficients, inequation sign and constant term of equation no %d seperated by spaces:\n" , i + 1);
        for(j = 0 ; j < in_var ; j++)
        {
            scanf("%f",&mat[i][j]);
        }
    cin >> inequality;
    if (inequality.at(0) == '<')
    {
        //cout << "here" << " " << i << " " << j << endl;
        mat[i][var] = 1;
        for(k = 0; k < eqn, k!=i ; k++)
        {
            mat[k][var] = 0;
        }
        var++;
    }
    else if (inequality.at(0) == '>')
    {
        mat[i][var] = 1;
        for(k = 0; k < eqn, k!=i ; k++)
        {
            mat[k][var] = 0;
        }
        var++;
        mat[i][var] = -1;
        for(k = 0; k < eqn, k!=i ; k++)
        {
            mat[k][var] = 0;
        }
        var++;
    }
    cin >> constants[i];
    }
    cout << "total variables formed: " << var << endl;
    for(i = 0 ; i < eqn ; i++)
    {
       mat[i][var] = constants[i];
    }
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

    /*for(i = 0 ; i < eqn ; i++)
    {
        for(j = 0 ; j <= var ; j++){
            cout << temp[i][j] << " ";
        }
    cout << endl;
    }*/

    rank2 = rankOfMatrix(eqn,var + 1);
    printf("Rank of A matrix is : %d \n", rank1);
    printf("Rank of A/b matrix is : %d \n", rank2);
    if(rank1 != rank2){
        printf("Cant solve\n");
        return 0;
    }
    if((rank1 = rank2) && (rank1 < eqn))
    {
       cout << "Equations are not linearly independent !" << endl;
        reduce(rank1, var, eqn);
      // cout << "Please remove the extra equation" << endl;
    //return 0;
    }
    for(int k = 0 ; k < rank1 ; k++)
    {
        for(j = 0 ; j <= var ; j++){
            cout << mat[k][j] << " ";
        }
        cout << endl;
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
        for(int k = 0 ; k < rank1 ; k++)
        {
            v = i;
            ind = 0;
            for(j = 0 ; j < var ; j++){
                if(v%2 == 1){
                    temp[k][ind++] = mat[k][j];
                }
                v = v/2;
            }
            temp[k][ind] = mat[k][var];
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
                }
                v = v/2;
            }
            temp[k][ind] = mat[k][var];
        }
        ind = 0;
        solution(var - m);
        v = i;
        
        for(j = 0 ; j < var ; j++){
            if(v%2 == 1){
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
        {
    countVal1++;
    cout<<"Non basic feasible solution is:";
    for(int l = 0 ; l < var ; l++){
            printf("%6.2f ", ans[k][l]);
        }
    cout << endl;
    continue;   
    }
        cout << countVal+1 << "th Basic feasible solution is:";
    countVal++;
    countVal1++;
        for(int l = 0 ; l < var ; l++){
            printf("%6.2f ", ans[k][l]);
        }
        printf("\n");
        
    }
    cout << "Total no of Basic Solutions are: " << countVal1 << endl;
    cout << "Total no of Basic Feasible Solutions are: " << countVal << endl;
    cout << "How many objective functions you want to optimize with these conditions" << endl;
    int parts=0;
    cin >> parts;

    for(int loop=0;loop<parts;loop++)
    {
        printf("Enter coefficients of each variable as in objective function\n");
        for(i = 0 ; i < in_var ; i++){
            scanf("%f",&z[i]);
        }
        printf("Enter 1 to get maximum, 2 to get minimum, 3 to get both \n");
        int bl = 0;
        scanf("%d", &bl);
        if(bl&(1 << 1))
            getMinimum(solInd, var);
        if(bl&1)
            getMaximum(solInd, var);
    }
    return 0;
}


