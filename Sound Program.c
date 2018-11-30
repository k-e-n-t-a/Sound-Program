// Sound_Program
// 任意の音楽データを作成する

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include "wave.h"

#define PI M_PI

#define ao_v(A)   ( 0.0 					   		 	) // a_0の式を代入
#define an_v(A,k) ( 4.0*A/(k*k*PI*PI)*(1.0-cos(k*PI))	) // a_nの式を代入
#define bn_v(A,k) ( 0.0								    ) // b_nの式を代入

#define ao_b(A)   ( 0.0 					   	 		) // a_0の式を代入
#define an_b(A,k) ( 0.0							 		) // a_nの式を代入
#define bn_b(A,k) ( -A/(k*PI)					 		) // b_nの式を代入

int main(void)
{
	MONO_PCM pcm1 ;
	int count_v, n_v, v[2][72] ;
	int count_b, n_b, b[7][27] ;
	int key, tempo, note, d, sum, n=0, n_j, j, p, c, k ;
	double tv=0.0, tb=0.0, y_b[7] ;
	double Amp, decayrate, saturation, A, ts, fs, fo, f[100] ;
	double y, Ao_v, Ao_b, An=0.0, Bn=0.0 ;
	
	v[0][0]=4;   v[1][0]=0;
	v[0][1]=8;   v[1][1]=76;
	v[0][2]=8;   v[1][2]=76;
	v[0][3]=8;   v[1][3]=74;
	v[0][4]=4;   v[1][4]=73;
	v[0][5]=4;   v[1][5]=74;
	v[0][6]=8;   v[1][6]=74;
	v[0][7]=8;   v[1][7]=73;
	v[0][8]=8;   v[1][8]=73;
	v[0][9]=8;   v[1][9]=69;
	v[0][10]=8;  v[1][10]=64;
	v[0][11]=4;  v[1][11]=64;
	
	v[0][12]=4;  v[1][12]=0;
	v[0][13]=8;  v[1][13]=76;
	v[0][14]=8;  v[1][14]=76;
	v[0][15]=8;  v[1][15]=74;
	v[0][16]=4;  v[1][16]=73;
	v[0][17]=4;  v[1][17]=74;
	v[0][18]=4;  v[1][18]=74;
	v[0][19]=8;  v[1][19]=73;
	v[0][20]=8;  v[1][20]=70;
	v[0][21]=8;  v[1][21]=71;
	v[0][22]=8;  v[1][22]=73;
	v[0][23]=8;  v[1][23]=71;
	
	v[0][24]=8;  v[1][24]=0;
	v[0][25]=8;  v[1][25]=74;
	v[0][26]=8;  v[1][26]=74;
	v[0][27]=8;  v[1][27]=74;
	v[0][28]=8;  v[1][28]=76;
	v[0][29]=8;  v[1][29]=74;
	v[0][30]=8;  v[1][30]=73;
	v[0][31]=8;  v[1][31]=74;
	v[0][32]=8;  v[1][32]=0;
	v[0][33]=8;  v[1][33]=74;
	v[0][34]=8;  v[1][34]=74;
	v[0][35]=8;  v[1][35]=74;
	v[0][36]=8;  v[1][36]=76;
	v[0][37]=8;  v[1][37]=74;
	v[0][38]=8;  v[1][38]=73;
	v[0][39]=8;  v[1][39]=74;
	
	v[0][40]=8;  v[1][40]=0;
	v[0][41]=8;  v[1][41]=71;
	v[0][42]=8;  v[1][42]=71;
	v[0][43]=8;  v[1][43]=71;
	v[0][44]=8;  v[1][44]=71;
	v[0][45]=8;  v[1][45]=71;
	v[0][46]=8;  v[1][46]=73;
	v[0][47]=8;  v[1][47]=74;
	v[0][48]=8;  v[1][48]=78;
	v[0][49]=8;  v[1][49]=76;
	v[0][50]=8;  v[1][50]=76;
	v[0][51]=8;  v[1][51]=73;
	v[0][52]=8;  v[1][52]=76;
	v[0][53]=8;  v[1][53]=76;
	v[0][54]=8;  v[1][54]=78;
	v[0][55]=8;  v[1][55]=76;
	
	
	b[0][0]=1; b[1][0]=67; b[2][0]=59; b[3][0]=55; b[4][0]=50; b[5][0]=47; b[6][0]=43;
	b[0][1]=1; b[1][1]=66; b[2][1]=59; b[3][1]=55; b[4][1]=50; b[5][1]=47; b[6][1]=43;
	b[0][2]=1; b[1][2]=65; b[2][2]=60; b[3][2]=57; b[4][2]=53; b[5][2]=48; b[6][2]=41;
	b[0][3]=1; b[1][3]=64; b[2][3]=59; b[3][3]=56; b[4][3]=52; b[5][3]=47; b[6][3]=40;
	
	b[0][4]=1; b[1][4]=64; b[2][4]=60; b[3][4]=57; b[4][4]=52; b[5][4]=45; b[6][4]=0;
	b[0][5]=1; b[1][5]=64; b[2][5]=60; b[3][5]=56; b[4][5]=52; b[5][5]=45; b[6][5]=0;
	b[0][6]=1; b[1][6]=64; b[2][6]=60; b[3][6]=55; b[4][6]=52; b[5][6]=45; b[6][6]=0;
	b[0][7]=1; b[1][7]=66; b[2][7]=62; b[3][7]=57; b[4][7]=50; b[5][7]=0;  b[6][7]=0;
	
	
	count_v = 55;									// v[][j]のjを設定
	count_b = 7;									// b[][j]のjを設定
	
	tempo = 68;										// テンポの設定（1分間の4分音符の数）
	key = 3;										// ここでキーの調節
	Amp = 1.0;										// 音の大きさを設定
	decayrate = 1.0001;								// ここに音の減衰率を設定
	saturation = 0.2;								// ここに減衰した音の最小値を設定
	
	for(j= 0; j<= count_v; j++){ tv=tv+1.0/(v[0][j]); }
	n_v=4.0*44100.0*60.0/(double)tempo*(double)tv;  // vocal
	
	for(j= 0; j<= count_b; j++){ tb=tb+1.0/(b[0][j]); }
	n_b=4.0*44100.0*60.0/(double)tempo*(double)tb;  // base
	
	if(n_v > n_b)  sum=n_v;
	if(n_b >= n_v) sum=n_b;
	
	printf("\n\t\tデータ数n[個]\t時間t[s]\n");
	printf("ヴォーカル:\t%d\t\t%d\n",n_v,n_v/44100);
	printf("ベース:\t\t%d\t\t%d\n",n_b,n_b/44100);
	printf("\n音楽ファイル（wav形式）を作成します。\n");
	getch();
	printf("\n\n\t\t…計算中です…\n");
	
	pcm1.fs = 44100; 								/* 標本化周波数 */
    pcm1.bits = 16; 								/* 量子化精度 */
	pcm1.length = sum+10; 							/* 音データの長さ */
	pcm1.s = calloc(pcm1.length, sizeof(double)); 	/* 音データ */
	
	ts=1.0/pcm1.fs;
	Ao_v=ao_v(A);
	Ao_b=ao_b(A);
	
	f[0]=0.0;										// 音の高さの初期設定始め
	for(d=20;d<=90;d++){
		f[d]=pow(2.0,((double)d-69.0+key)/12.0)*440.0;
	}
	
	for(j=0;j<=count_v;j++){
		
		note=v[0][j];
		d=v[1][j];
		fo=f[d];
		n_j=4.0*pcm1.fs/(double)note*60.0/(double)tempo;
		
		//printf("\n\nヴォーカル：[%d]\t音符\tn[個]\tt[s]\t音階d\n",j);
		//printf("\t\t%d\t%d\t%.3lf\t%d\n",note,n_j,n_j*ts,d);
		//printf("\nn[個]\tt[s]\ty\td\tfo[Hz]\t振幅A\tpcm1.s[n]\n");
		//getch();
		
		for(p=0;p<=n_j;p++){						/* ts刻みで音データを格納 */
			
			A=Amp*pow(decayrate,-(double)p)+saturation;
			
			for(k= 1 ;k<= 10 ;k++){
				
				An+=an_v(A,k)*cos((double)k*2.0*PI*fo*n*ts);
				Bn+=bn_v(A,k)*sin((double)k*2.0*PI*fo*n*ts);
			}
			
			y=Ao_v+An+Bn;
			pcm1.s[n]=y;								/* 音データのコピー */
			//printf("%d\t%.3lf\t%.3lf\t%d\t%.1lf\t%.3lf\t%.3lf\n",n,n*ts,y,d,fo,A,pcm1.s[n]);
			n=n+1;
			An=0.0;
			Bn=0.0;
		}
		//printf("n[個]\tt[s]\ty\td\tfo[Hz]\t振幅A\tpcm1.s[n]\n");
	}
	
	n=0;
	
	for(j=0;j<=count_b;j++){
		
		note=b[0][j];
		n_j=4.0*pcm1.fs/(double)note*60.0/(double)tempo;
		
		
		//printf("\n\nベース：[%d]\t音符\tn[個]\tt[s]\t音階d\n",j);
		//printf("\t\t%d\t%d\t%.3lf\t%d,%d,%d,%d,%d,%d\n",note,n_j,n_j*ts,b[6][j],b[5][j],b[4][j],b[3][j],b[2][j],b[1][j]);
		//printf("\nn[個]\tt[s]\ty\t振幅A\tpcm1.s[n]\n");
		//getch();
		
		
		for(p=0;p<=n_j;p++){						/* ts刻みで音データを格納 */
			
			A=Amp*pow(decayrate,-(double)p)+saturation;
			
			for(c=1;c<=6;c++){
				
				d=b[c][j];
				fo=f[d];
				
				for(k= 1 ;k<= 10 ;k++){
				
					An+=an_b(A,k)*cos((double)k*2.0*PI*fo*n*ts);
					Bn+=bn_b(A,k)*sin((double)k*2.0*PI*fo*n*ts);
				}
				
				y_b[c]=Ao_b+An+Bn;
				An=0.0;
				Bn=0.0;
			}
			y=y_b[1]+y_b[2]+y_b[3]+y_b[4]+y_b[5]+y_b[6];  /* 音データのコピー */
			pcm1.s[n]=pcm1.s[n]+y;
			
			//printf("%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",n,n*ts,y,A,pcm1.s[n]);
			n=n+1;
		}
		//printf("n[個]\tt[s]\ty\t振幅A\tpcm1.s[n]\n");
	}
	
	wave_write_16bit_mono(&pcm1, "Sound_Program.wav");		/* 音データの出力 */
    free(pcm1.s);
	
	printf("\n\nファイル「Sound_Program.wav」に保存しました。\n");
	printf("\nデータ数n:\t%d\t[個]\n時間t:\t\t%d\t[s]\n\n",pcm1.length,pcm1.length/pcm1.fs);
	
	return(0);
}
