

void sparse_set(double *Pr,int *Ir,int *Jc,int ri,int rj,double val)
{
// Binary search within column

	int k,k1,k2,cr=ri;
	k1=Jc[rj];
	k2=Jc[rj+1]-1;
	do{
		k=(k1+k2)>>1;
		if(cr<Ir[k])
			k2=k-1;
		else
			k1=k+1;
		if(cr==Ir[k])
			break;
	}while(k2>=k1);

	/*if(ri!=Ir[k]){
		char str[50];
			sprintf(str,"fprintf\('%d,%d,'\)",ri,rj);
			mexEvalString(str);
			return;
	}*/
	Pr[k]+=val;
}
