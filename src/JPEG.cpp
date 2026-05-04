#include "JPEG.h"
#include "NxNDCT.h"
#include <math.h>

#include "JPEGBitStreamWriter.h"


#define DEBUG(x) do{ qDebug() << #x << " = " << x;}while(0)



// quantization tables from JPEG Standard, Annex K
uint8_t QuantLuminance[8*8] =
    { 16, 11, 10, 16, 24, 40, 51, 61,
      12, 12, 14, 19, 26, 58, 60, 55,
      14, 13, 16, 24, 40, 57, 69, 56,
      14, 17, 22, 29, 51, 87, 80, 62,
      18, 22, 37, 56, 68,109,103, 77,
      24, 35, 55, 64, 81,104,113, 92,
      49, 64, 78, 87,103,121,120,101,
      72, 92, 95, 98,112,100,103, 99 };
uint8_t QuantChrominance[8*8] =
    { 17, 18, 24, 47, 99, 99, 99, 99,
      18, 21, 26, 66, 99, 99, 99, 99,
      24, 26, 56, 99, 99, 99, 99, 99,
      47, 66, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99 };

static char quantizationMatrix[64] =
{
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

struct imageProperties{
    int width;
    int height;
    int16_t* coeffs;
};


void DCTUandV(const char input[], int16_t output[], int N, double* DCTKernel)
{
    double* temp = new double[N*N];
    double* DCTCoefficients = new double[N*N];

    double sum;
    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = 0; j <= N - 1; j++)
        {
            sum = 0;
            for (int k = 0; k <= N - 1; k++)
            {
                sum = sum + DCTKernel[i*N+k] * (input[k*N+j]);
            }
            temp[i*N + j] = sum;
        }
    }

    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = 0; j <= N - 1; j++)
        {
            sum = 0;
            for (int k = 0; k <= N - 1; k++)
            {
                sum = sum + temp[i*N+k] * DCTKernel[j*N+k];
            }
            DCTCoefficients[i*N+j] = sum;
        }
    }

    for(int i = 0; i < N*N; i++)
    {
        output[i] = floor(DCTCoefficients[i]+0.5);
    }

    delete[] temp;
    delete[] DCTCoefficients;

    return;
}

uint8_t quantQuality(uint8_t quant, uint8_t quality) {
    // Convert to an internal JPEG quality factor, formula taken from libjpeg
    int16_t q = quality < 50 ? 5000 / quality : 200 - quality * 2;
    return clamp((quant * q + 50) / 100, 1, 255);
}

static void doZigZag(int16_t block[], uint8_t quantizationBlock[], int N, int DCTorQuantization)
{
    if(DCTorQuantization)  // DCT
        {
            int16_t* temp = new int16_t[N*N];
            int x = 0, y = 0;
            bool up = true;
            for (int i = 0; i < N*N; i++)
            {
                temp[i] = block[y*N + x];
                if(up){
                    if(x == N-1)      { y++;  up = false; }     // udarimo desno
                    else if(y == 0)   { x++;  up = false; }     // udarimo gore
                    else              { x++;  y--;        }     // krecemo se dijagonalno gore
                } else {
                    if(y == N-1)      { x++;  up = true;  }     // udarimo dole
                    else if(x == 0)   { y++;  up = true;  }     // udarimo levo
                    else              { x--;  y++;        }     // krecemo se dijagonalno dole
                }
            }
            for(int i = 0; i < N*N; i++) block[i] = temp[i];
            delete[] temp;
        }
    else  // Quantization
        {
            uint8_t* temp = new uint8_t[N*N];
            int x = 0, y = 0;
            bool up = true;
            for (int i = 0; i < N*N; i++)
            {
                temp[i] = quantizationBlock[y*N + x];
                if(up){
                    if(x == N-1)      { y++;  up = false; }
                    else if(y == 0)   { x++;  up = false; }
                    else              { x++;  y--;        }
                } else {
                    if(y == N-1)      { x++;  up = true;  }
                    else if(x == 0)   { y++;  up = true;  }
                    else              { x--;  y++;        }
                }
            }
            for(int i = 0; i < N*N; i++) quantizationBlock[i] = temp[i];
            delete[] temp;
        }
}

/* perform DCT */
imageProperties performDCT(char input[], int xSize, int ySize, int N, uint8_t quality, bool quantType)
{
    imageProperties output;

    output.width = xSize/ N;
    output.height = ySize / N;

    output.coeffs = new int16_t[xSize * ySize];

    double *kernel = new double[N*N];
    GenerateDCTmatrix(kernel,N);

    uint8_t quantTable[64];
    uint8_t *baseTable = quantType ? QuantLuminance : QuantChrominance;
    for(int i = 0; i < N*N; i++){
        quantTable[i] = quantQuality(baseTable[i],quality);
    }

    for(int by = 0; by < output.height; by++){
        for(int bx = 0; bx < output.width; bx++){
            char block[64];

            for(int row = 0; row < N; row++)
                for(int col = 0; col < N; col++)
                    block[row*N + col] = input[(by*N + row)*xSize + (bx*N + col)];  // direktno pristupamo bloku (br.bloka(vertikalnog ili horizontalnog) * vel.bloka) + odgovarajuci red ili kolona

            int16_t *blockCoeffs = &output.coeffs[(by*output.width + bx) * N * N];  // pokazivac na lokaciju gdje treba da upisemo koeficijente u memoriji (mnozimo sa 64 jer nam je jedan blok toliki)
            DCTUandV(block, blockCoeffs, N, kernel);

            for(int i = 0; i < N*N; i++)
                blockCoeffs[i] = round(blockCoeffs[i] / quantTable[i]);

            doZigZag(blockCoeffs,nullptr,N,1);
        }
    }

    delete[] kernel;
    return output;
}

//JPEGBitStreamWriter streamer("example.jpg");
void performJPEGEncoding(uchar Y_buff[], char U_buff[], char V_buff[], int xSize, int ySize, int quality)
{
	DEBUG(quality);
	
	
    auto s = new JPEGBitStreamWriter("example.jpg");
	// TO DO
    const int N = 8;

    // pomeramo opseg za Y komponentu
    char *Y_temp = new char[xSize * ySize];
    for(int i = 0; i < xSize * ySize; i++){
        Y_temp[i] = char(Y_buff[i] - 128);
    }

    // vrsimo prosirivanje
    char *Y_extended, *U_extended, *V_extended;
    int newXSize, newYSize, notNeeded;

    extendBorders(Y_temp, xSize, ySize, N*2, &Y_extended, &newXSize, &newYSize);
    extendBorders(U_buff, xSize/2, ySize/2, N, &U_extended, &notNeeded, &notNeeded);
    extendBorders(V_buff, xSize/2, ySize/2, N, &V_extended, &notNeeded, &notNeeded);

    delete[] Y_temp;

    imageProperties Y_properties = performDCT(Y_extended,newXSize,newYSize,N,quality,1);
    imageProperties U_properties = performDCT(U_extended,newXSize/2,newYSize/2,N,quality,0);
    imageProperties V_properties = performDCT(V_extended,newXSize/2,newYSize/2,N,quality,0);

    // kvantizacione matrice
    uint8_t lumiQuant[64], chromaQuant[64];

    for(int i = 0; i < N*N; i++){
        lumiQuant[i] = quantQuality(QuantLuminance[i],quality);
        chromaQuant[i] = quantQuality(QuantChrominance[i],quality);
    }

    uint8_t *zigZagLumi = new uint8_t[N*N];
    uint8_t *zigZagChroma = new uint8_t[N*N];

    for(int i = 0; i < N*N; i++){
        zigZagLumi[i] = lumiQuant[i];
        zigZagChroma[i] = chromaQuant[i];
    }

    // zig zag kvantizacionih matrica
    doZigZag(nullptr,zigZagLumi,N,0);
    doZigZag(nullptr,zigZagChroma,N,0);

    s->writeHeader();
    s->writeQuantizationTables(zigZagLumi,zigZagChroma);
    s->writeImageInfo(xSize,ySize);
    s->writeHuffmanTables();

    // broj 16*16 blokova
    int macroBlocksX = newXSize / 16;
    int macroBlocksY = newYSize / 16;

    for(int mby = 0; mby < macroBlocksY; mby++)
    {
        for(int mbx = 0; mbx < macroBlocksX; mbx++)
        {
            int yW = Y_properties.width;  // broj blokova u redu za Y
            s->writeBlockY(&Y_properties.coeffs[( mby*2    * yW + mbx*2   ) * 64]);
            s->writeBlockY(&Y_properties.coeffs[( mby*2    * yW + mbx*2+1 ) * 64]);
            s->writeBlockY(&Y_properties.coeffs[((mby*2+1) * yW + mbx*2   ) * 64]);
            s->writeBlockY(&Y_properties.coeffs[((mby*2+1) * yW + mbx*2+1 ) * 64]);

            int uvW = U_properties.width;  // ... za U i V
            s->writeBlockU(&U_properties.coeffs[(mby * uvW + mbx) * 64]);

            s->writeBlockV(&V_properties.coeffs[(mby * uvW + mbx) * 64]);
        }
    }

    s->finishStream();

    delete[] zigZagChroma;
    delete[] zigZagLumi;

}
