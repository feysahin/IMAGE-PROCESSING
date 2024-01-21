#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>


/*
		Görseller, "IMAGES" klasöründe bulunan kendilerine ait klasörlere kaydedilmektedir. (IMAGES> lena|fruit|coins)
*/

typedef struct pgmData {
    int row;
    int col;
    int max_color;
    char pgmType[2];
    unsigned char **matrix;
}PGM;

void printMatrix(unsigned char **m, int r, int c);
unsigned char **matrixAllocate(int r, int c);
int **matrixAllocate_int(int r, int c);
void *get_img(char fileName[], PGM *pgm);
void convolution(int **after_convolution, unsigned char **pgmMatrix, int r, int c, int k_size, float kernel[][k_size]);
void create_gauss_kernel(int size, float gaussKernel[][size], float sigma);
void skipComments(FILE* fp);
void skipSpace(FILE* fp);
void save_file(char fileName[], PGM* pgm, int **after_convolution);
void minmax_normalization(int **after_convolution, int r, int c);
void sobel_filter(int **after_xy, int **after_x, int **after_y, PGM *pgm);
void laplacian_filter(int **after_convolution, PGM *pgm, float lap[][3]);

int main(){
	
	int size = 7;
	float sigma = 4;
	
	char read_fName[] = "lena.pgm", 
		 // gauss
		 gauss_fName[] = "IMAGES/lena/lena_gauss.pgm", 
		 // sobel without gauss
		 x_sobel_fName[] = "IMAGES/lena/lena_sobel_x.pgm", 
		 y_sobel_fName[] = "IMAGES/lena/lena_sobel_y.pgm", 
		 xy_sobel_fName[] = "IMAGES/lena/lena_sobel_xy.pgm", 
		 // sobel after gauss
		 gauss_x_sobel_fName[] = "IMAGES/lena/lena_gauss_sobel_x.pgm", 
		 gauss_y_sobel_fName[] = "IMAGES/lena/lena_gauss_sobel_y.pgm", 
		 gauss_xy_sobel_fName[] = "IMAGES/lena/lena_gauss_sobel_xy.pgm",
		 // laplacian after gauss
		 gauss_lap1_fName[] = "IMAGES/lena/lena_gauss_lap1.pgm",
		 gauss_lap2_fName[] = "IMAGES/lena/lena_gauss_lap2.pgm";
		 
	
	/*char read_fName[] = "fruit.pgm", 
		 gauss_fName[] = "IMAGES/fruit/fruit_gauss.pgm", 
		 x_sobel_fName[] = "IMAGES/fruit/fruit_sobel_x.pgm", 
		 y_sobel_fName[] = "IMAGES/fruit/fruit_sobel_y.pgm", 
		 xy_sobel_fName[] = "IMAGES/fruit/fruit_sobel_xy.pgm", 
		 gauss_x_sobel_fName[] = "IMAGES/fruit/fruit_gauss_sobel_x.pgm", 
		 gauss_y_sobel_fName[] = "IMAGES/fruit/fruit_gauss_sobel_y.pgm", 
		 gauss_xy_sobel_fName[] = "IMAGES/fruit/fruit_gauss_sobel_xy.pgm",
		 gauss_lap1_fName[] = "IMAGES/fruit/fruit_gauss_lap1.pgm",
		 gauss_lap2_fName[] = "IMAGES/fruit/fruit_gauss_lap2.pgm";*/
		 
	/*
	char read_fName[] = "coins_ascii.pgm", 
		 gauss_fName[] = "IMAGES/coins/coins_gauss.pgm", 
		 x_sobel_fName[] = "IMAGES/coins/coins_sobel_x.pgm", 
		 y_sobel_fName[] = "IMAGES/coins/coins_sobel_y.pgm", 
		 xy_sobel_fName[] = "IMAGES/coins/coins_sobel_xy.pgm", 
		 gauss_x_sobel_fName[] = "IMAGES/coins/coins_gauss_sobel_x.pgm", 
		 gauss_y_sobel_fName[] = "IMAGES/coins/coins_gauss_sobel_y.pgm", 
		 gauss_xy_sobel_fName[] = "IMAGES/coins/coins_gauss_sobel_xy.pgm",
		 gauss_lap1_fName[] = "IMAGES/coins/coins_gauss_lap1.pgm",
		 gauss_lap2_fName[] = "IMAGES/coins/coins_gauss_lap2.pgm";*/
		 	
	// READ FILE
	PGM pgm;
	get_img(read_fName, &pgm);
	
	
	// SOBEL 
	int **after_x, **after_y, **after_xy;
	after_x = matrixAllocate_int(pgm.row, pgm.col);
	after_y = matrixAllocate_int(pgm.row, pgm.col);
	after_xy = matrixAllocate_int(pgm.row, pgm.col);
	
	// before gauss
	/*sobel_filter(after_xy, after_x, after_y, &pgm);
	save_file(x_sobel_fName, &pgm, after_x);
	save_file(y_sobel_fName, &pgm, after_y);
	save_file(xy_sobel_fName, &pgm, after_xy);*/

	
	// GAUSSIAN FILTER
	float gaussKernel[size][size];
	create_gauss_kernel(size, gaussKernel, sigma);
	
	int **after_gauss;
	after_gauss = matrixAllocate_int(pgm.row, pgm.col);
	convolution(after_gauss, pgm.matrix, pgm.row, pgm.col, size, gaussKernel);
	minmax_normalization(after_gauss, pgm.row, pgm.col);
	save_file(gauss_fName, &pgm, after_gauss);
	
		
	// after gauss
	int i, j;
	for(i=0; i<pgm.row; i++){
		for(j=0; j<pgm.col; j++){
			after_x[i][j] = (int) after_gauss[i][j];
			after_y[i][j] = (int) after_gauss[i][j];
			after_xy[i][j] = (int) after_gauss[i][j];
		}
	}
	sobel_filter(after_xy, after_x, after_y, &pgm);
	save_file(gauss_x_sobel_fName, &pgm, after_x);
	save_file(gauss_y_sobel_fName, &pgm, after_y);
	save_file(gauss_xy_sobel_fName, &pgm, after_xy);
	
	
	// LAPLACIAN
	int **after_lap1, **after_lap2;
	after_lap1 = matrixAllocate_int(pgm.row, pgm.col);
	after_lap2 = matrixAllocate_int(pgm.row, pgm.col);

	for(i=0; i<pgm.row; i++){
		for(j=0; j<pgm.col; j++){
			after_lap1[i][j] = (int) after_gauss[i][j];
			after_lap2[i][j] = (int) after_gauss[i][j];
		}
	}
	// laplacian-1
	float lap1[3][3] = {{0, -1, 0}, {-1, 4, -1}, {0, -1, 0}};
	laplacian_filter(after_lap1, &pgm, lap1);
	save_file(gauss_lap1_fName, &pgm, after_lap1);
	
	// laplacian-2
	float lap2[3][3] = {{-1, -1, -1}, {-1, 8, -1}, {-1, -1, -1}};
	laplacian_filter(after_lap2, &pgm, lap2);
	save_file(gauss_lap2_fName, &pgm, after_lap2);
	
	
	free(pgm.matrix);
	free(after_gauss);
	
	free(after_x);
	free(after_y);
	free(after_xy);
	
	free(after_lap1);
	free(after_lap2);
	
	return 0;
}

void laplacian_filter(int **after_convolution, PGM *pgm, float lap[][3]){
	
	int **after_lap, min = 1000, max = -1000, i, j, k, p, sum;
	after_lap = matrixAllocate_int(pgm->row, pgm->col);
		
	// Convolution
	for(i = 0; i < pgm->row - 2; i++){
		for(j = 0; j < pgm->col - 2; j++){
			sum = 0;
			for(k = 0; k < 3; k++){
				for(p = 0; p < 3; p++){
					sum += after_convolution[i+k][j+p] * lap[k][p];
				}
			}
			after_lap[i + 1][j + 1] = sum;
			
			if(sum > max) max = sum;
            else if(sum < min) min = sum;
		}
	}
	// Normalization
	for (i = 0; i < pgm->row; i++){
        for (j = 0; j < pgm->col; j++){
        	after_lap[i][j] = (int) (((float) (after_lap[i][j] - min) / (max - min)) * 255);
        }
    }
	// Copy
	for(i=0; i<pgm->row; i++){
		for(j=0; j<pgm->col; j++){
			after_convolution[i][j] = after_lap[i][j];
		}
	}
}

void sobel_filter(int **after_xy, int **after_x, int **after_y, PGM *pgm){
	
	int i, j, k, p;
	float sobel_x[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}; 
	float sobel_y[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};	
						 
	int **x_sobeled, **y_sobeled;
	x_sobeled = matrixAllocate_int(pgm->row, pgm->col);
	y_sobeled = matrixAllocate_int(pgm->row, pgm->col);

	// copy of original for padding
	for(i=0; i<pgm->row; i++){
		for(j=0; j<pgm->col; j++){
			x_sobeled[i][j] = (int) pgm->matrix[i][j];
			y_sobeled[i][j] = (int) pgm->matrix[i][j];
		}
	}
	int x_min = 1000, x_max = -1000, x_sum;
	int y_min = 1000, y_max = -1000, y_sum;
	
	// Convolution
	for(i = 0; i < pgm->row - 2; i++){
		for(j = 0; j < pgm->col - 2; j++){
			x_sum = 0;
			y_sum = 0;
			for(k = 0; k < 3; k++){
				for(p = 0; p < 3; p++){
					x_sum += pgm->matrix[i+k][j+p] * sobel_x[k][p];
					y_sum += pgm->matrix[i+k][j+p] * sobel_y[k][p];
				}
			}
			x_sobeled[i + 1][j + 1] = x_sum;
			y_sobeled[i + 1][j + 1] = y_sum;
			
			if (x_sum > x_max) x_max = x_sum;
            else if (x_sum < x_min)	x_min = x_sum;
            
            if (y_sum > y_max) y_max = y_sum;
            else if (y_sum < y_min) y_min = y_sum;
		}
	}	
	// XY-sobel
    for(i = 0; i < pgm->row; i++){
		for(j = 0; j < pgm->col; j++){
			after_xy[i][j] = (int) (sqrt(pow(x_sobeled[i][j], 2) + pow(y_sobeled[i][j], 2)));
		}
	}
	// normalization
	for (i = 0; i < pgm->row; i++){
        for (j = 0; j < pgm->col; j++){
        	x_sobeled[i][j] =  (int) (((float) (x_sobeled[i][j] - x_min) / (x_max - x_min)) * 255);
        	y_sobeled[i][j] =  (int) (((float) (y_sobeled[i][j] - y_min) / (y_max - y_min)) * 255);
        }
    }
	// X and Y sobels
	for(i=0; i<pgm->row; i++){
		for(j=0; j<pgm->col; j++){
			after_x[i][j] = x_sobeled[i][j];
			after_y[i][j] = y_sobeled[i][j];
		}
	}
	int min = 1000, max = -1000;
	for(i=0; i<pgm->row; i++){
		for(j=0; j<pgm->col; j++){
			if(after_xy[i][j] < min){
				min = after_xy[i][j];
			}
			else if(after_xy[i][j] > max){
				max = after_xy[i][j];
			}
		}
	}
	// normalization
	for(i=0; i<pgm->row; i++){
		for(j=0; j<pgm->col; j++){
			after_xy[i][j] =  (int) (((float) (after_xy[i][j] - min) / (max - min)) * 255);
		}
	}
	free(x_sobeled);
	free(y_sobeled);
}

void minmax_normalization(int **after_convolution, int r, int c){
	
	int i, j, min = after_convolution[0][0], max = after_convolution[0][0];
	
	for(i=0; i<r; i++){
		for(j=0; j<c; j++){
			if(after_convolution[i][j] < min){
				min = after_convolution[i][j];
			}
			else if(after_convolution[i][j] > max){
				max = after_convolution[i][j];
			}
		}
	}
	int dif = max-min;
	for(i=0; i<r; i++){
		for(j=0; j<c; j++){
			after_convolution[i][j] =  (int) (((float) (after_convolution[i][j] - min) / dif) * (255-min)) + min;
		}
	}
}

void convolution(int **after_convolution, unsigned char **pgmMatrix, int r, int c, int k_size, float kernel[][k_size]){
	
	int i, j, k, p, sum;
	
	// copy of Original for Padding
	for(i=0; i<r; i++){
		for(j=0; j<c; j++){
			after_convolution[i][j] = (int) pgmMatrix[i][j];
		}
	}
	// Convolution
	for(i = 0; i < r - k_size+1; i++){
		for(j = 0; j < c - k_size+1; j++){
			sum = 0;
			for(k = 0; k < k_size; k++){
				for(p = 0; p < k_size; p++){
					sum += pgmMatrix[i+k][j+p] * kernel[k][p];
				}
			}
			after_convolution[i + k_size/2][j + k_size/2] = round( (float) sum / (k_size*k_size));
		}
	}	
}

void *get_img(char fileName[], PGM *pgm){

	// P2-> fscanf 
	// P5-> fread
		
	FILE *file = fopen(fileName, "rb");
	if(file != NULL){
		int i, j;
	
		// P5/P2
		fscanf(file, "%s", pgm->pgmType); 
		printf("%s\n", pgm->pgmType);		
				
		// comments
		skipComments(file); 
		
		// row-col
		fscanf(file, "%d %d", &pgm->col, &pgm->row); 
		printf("row: %d\ncol: %d\n", pgm->row, pgm->col);
		
		// max_color
		fscanf(file, "%d", &pgm->max_color);
		printf("max color range: %d\n\n", pgm->max_color);
		
		pgm->matrix = matrixAllocate(pgm->row, pgm->col);
		
 		skipSpace(file);
 		// P2:
 		if(!strcmp(pgm->pgmType, "P2")){
 			for(i=0; i<pgm->row; i++){
 				for(j=0; j<pgm->col; j++){
	 				skipSpace(file);
	 				fscanf(file, "%d", &pgm->matrix[i][j]);
				}
		    }
		}
 		// P5:
		else{
			for(i=0; i<pgm->row; i++){
				fread(pgm->matrix[i], sizeof(unsigned char), pgm->col, file);
			}
		}
		fclose(file);
	}
	else printf("\nERROR! Could not open the file.\n");
}

void save_file(char fileName[], PGM* pgm, int **after_convolution){

	// P5: fwrite
	// P2: fprintf
	
	FILE* file = fopen(fileName, "wb");
	fprintf(file, "%s\n%d %d\n%d\n", pgm->pgmType, pgm->col, pgm->row, pgm->max_color);
	
	int i, j;
	for(i = 0; i < pgm->row; i++) {
		for(j = 0; j < pgm->col; j++){
			pgm->matrix[i][j] = (unsigned char) after_convolution[i][j];
		}
	}
	// P5
	if(!strcmp(pgm->pgmType, "P5")){
			
		for(i=0; i<pgm->row; i++){
			fwrite(pgm->matrix[i], sizeof(unsigned char), pgm->col, file);
		}
	} 
	// P2
	else{
		int i, j;
		for(i = 0; i < pgm->row; i++) {
			for(j = 0; j < pgm->col; j++){
				fprintf(file, "%d%c", (unsigned char) after_convolution[i][j], ((i+1) % pgm->col) ? ' ' : '\n');
			}
		}
	}
	fclose(file);
}

void create_gauss_kernel(int size, float gaussKernel[][size], float sigma){
	
	int i, j, X[size][size], Y[size][size], dif = (size-1)/2;
	float pi = 3.14, e = 2.73;
	
	// INITIALIZATION OF X-Y
	for(i=0; i<size; i++){		
		for(j=0; j<size; j++){
			X[i][j] = j-dif;
			Y[i][j] = i-dif;
		}
	}	
	// CREATE GAUSSIAN KERNEL
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			gaussKernel[i][j] = (pow(e, -( ( pow(X[i][j], 2) + pow(Y[i][j], 2)) / (2.0*pow(sigma, 2)))) / (2.00*pi*pow(sigma, 2)));
		}
	}
	float min = gaussKernel[0][0];
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			gaussKernel[i][j] = gaussKernel[i][j] / min;
		}
	}	
	printf("\n\n-------------------\nGaussian Kernel:\n\n");
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			printf("  %.3f ", gaussKernel[i][j]);
		}
		printf("\n");
	}	
}

void skipComments(FILE* fp){
    int ch;
    char line[100];
 
    while((ch = fgetc(fp)) != EOF && isspace(ch));

    if(ch == '#'){
        fgets(line, sizeof(line), fp);
        skipComments(fp);
    }
    else fseek(fp, -1, SEEK_CUR);
}

void skipSpace(FILE* fp){
    int ch; 
    while((ch = fgetc(fp)) != EOF && isspace(ch));
    if(ch != EOF) fseek(fp, -1, SEEK_CUR);
}

void printMatrix(unsigned char **m, int r, int c){
	int i, j;
	for(i=0; i<r; i++){
		for(j=0; j<c; j++){
			printf("%d ", m[i][j]);
		}
		printf("\n\n");
	}
}

unsigned char **matrixAllocate(int r, int c){
	int i;
	unsigned char **alloc_m;
	alloc_m = (unsigned char**)malloc(r*sizeof(unsigned char*));
	for(i=0; i<r; i++){
		alloc_m[i] = (unsigned char*)malloc(c*sizeof(unsigned char));	
	}
	return alloc_m;
}

int **matrixAllocate_int(int r, int c){
	int i;
	int **alloc_m;
	alloc_m = (int**)malloc(r*sizeof(int*));
	for(i=0; i<r; i++){
		alloc_m[i] = (int*)malloc(c*sizeof(int));	
	}
	return alloc_m;
}

