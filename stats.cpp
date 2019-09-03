
#include "stats.h"

stats::stats(PNG & im){

/* your code here */
    int height = im.height();
    int width = im.width();
    sumSat.resize(width, vector< double >(height));
    sumLum.resize(width, vector< double >(height));
    sumHueX.resize(width, vector< double >(height));
    sumHueY.resize(width, vector< double >(height));
    hist.resize(width, vector<vector<int>>(height, vector<int>(36)));

    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            HSLAPixel* pixel = im.getPixel(i, j);
            int k = (int)floor((pixel->h)/10);
            if (i == 0 && j == 0){
                sumSat[i][j] = pixel->s;
                sumLum[i][j] = pixel->l;
                sumHueX[i][j] = cos((pixel->h)*PI/180);
                sumHueY[i][j] = sin((pixel->h)*PI/180);
                hist[i][j][k]++;
            }else if (i == 0 && j > 0){
                sumSat[i][j] = sumSat[i][j-1] + pixel->s;
                sumLum[i][j] = sumLum[i][j-1] + pixel->l;
                sumHueX[i][j] = sumHueX[i][j-1] + cos((pixel->h)*PI/180);
                sumHueY[i][j] = sumHueY[i][j-1] + sin((pixel->h)*PI/180);
                for (int h = 0; h < 36; h++)
                    hist[i][j][h] = hist[i][j-1][h];
                hist[i][j][k] ++;
            }else if (i > 0 && j == 0){
                sumSat[i][j] = sumSat[i-1][j] + pixel->s;
                sumLum[i][j] = sumLum[i-1][j] + pixel->l;
                sumHueX[i][j] = sumHueX[i-1][j] + cos((pixel->h)*PI/180);
                sumHueY[i][j] = sumHueY[i-1][j] + sin((pixel->h)*PI/180);
                for (int h = 0; h < 36; h++)
                    hist[i][j][h] = hist[i-1][j][h];
                hist[i][j][k] ++;
            }else{
                sumSat[i][j] = sumSat[i-1][j] + sumSat[i][j-1] - sumSat[i-1][j-1] + pixel->s;
                sumLum[i][j] = sumLum[i-1][j] + sumLum[i][j-1] - sumLum[i-1][j-1] + pixel->l;
                sumHueX[i][j] = sumHueX[i-1][j] + sumHueX[i][j-1] - sumHueX[i-1][j-1] + cos((pixel->h)*PI/180);
                sumHueY[i][j] = sumHueY[i-1][j] + sumHueY[i][j-1] - sumHueY[i-1][j-1] + sin((pixel->h)*PI/180);
                for (int h = 0; h < 36; h++)
                    hist[i][j][h] = hist[i][j-1][h] + hist[i-1][j][h] - hist[i-1][j-1][h];
                hist[i][j][k] ++;
            }
        }
    }
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    return (lr.first - ul.first + 1)*(lr.second - ul.second + 1);

}

HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    long area = rectArea(ul, lr);
    HSLAPixel pixel = HSLAPixel();
    double averageHueX;
    double averageHueY;
    if (ul.first == 0 && ul.second == 0){
        pixel.s = sumSat[lr.first][lr.second]/area;
        pixel.l = sumLum[lr.first][lr.second]/area;
        averageHueX = sumHueX[lr.first][lr.second]/area;
        averageHueY = sumHueY[lr.first][lr.second]/area;
    }else if (ul.first == 0){
        pixel.s = (sumSat[lr.first][lr.second] - sumSat[lr.first][ul.second-1])/area;
        pixel.l = (sumLum[lr.first][lr.second] - sumLum[lr.first][ul.second-1])/area;
        averageHueX = (sumHueX[lr.first][lr.second] - sumHueX[lr.first][ul.second-1])/area;
        averageHueY = (sumHueY[lr.first][lr.second] - sumHueY[lr.first][ul.second-1])/area;
    }else if (ul.second == 0){
        pixel.s = (sumSat[lr.first][lr.second] - sumSat[ul.first-1][lr.second])/area;
        pixel.l = (sumLum[lr.first][lr.second] - sumLum[ul.first-1][lr.second])/area;
        averageHueX = (sumHueX[lr.first][lr.second] - sumHueX[ul.first-1][lr.second])/area;
        averageHueY = (sumHueY[lr.first][lr.second] - sumHueY[ul.first-1][lr.second])/area;
    }else{
        pixel.s = (sumSat[lr.first][lr.second] - sumSat[ul.first-1][lr.second] - sumSat[lr.first][ul.second-1] + sumSat[ul.first-1][ul.second-1 ])/area;
        pixel.l = (sumLum[lr.first][lr.second] - sumLum[ul.first-1][lr.second] - sumLum[lr.first][ul.second-1] + sumLum[ul.first-1][ul.second-1])/area;
        averageHueX = (sumHueX[lr.first][lr.second] - sumHueX[ul.first-1][lr.second] - sumHueX[lr.first][ul.second-1] + sumHueX[ul.first-1][ul.second-1])/area;
        averageHueY = (sumHueY[lr.first][lr.second] - sumHueY[ul.first-1][lr.second] - sumHueY[lr.first][ul.second-1] + sumHueY[ul.first-1][ul.second-1])/area;
    }
    double hue = atan2(averageHueY, averageHueX)* 180 / PI;
    if (hue < 0)
        hue += 360;
    pixel.h = hue;
    pixel.a = 1.0;
    return pixel;
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    vector<int> distn;
    for(int i = 0; i < 36; i++){
            distn.push_back(0);
        }
    if(lr.second < ul.second || lr.first<ul.first || lr.second<0 || lr.first<0 || ul.first<0 || ul.second<0){
        return distn;
    }
    for (int i = 0; i < 36; i++){
        distn[i] = hist[lr.first][lr.second][i];
    }
    if (ul.first != 0){
        for (int i = 0; i < 36; i++){
            distn[i] = distn[i] - hist[ul.first-1][lr.second][i];
        }
    }
    if(ul.second != 0){
        for (int i = 0; i < 36; i++){
            distn[i] = distn[i] - hist[lr.first][ul.second-1][i];
        }
        if (ul.first) {
			for (int i = 0; i < 36; i++)
				distn[i] = distn[i] + hist[ul.first - 1][ul.second - 1][i];
		}
    }
    return distn;
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){

    double entropy = 0.;

/* your code here */

    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 ) {
            entropy += ((double) distn[i]/(double) area) 
                                    * log2((double) distn[i]/(double) area);
        }
    }

    return  -1 * entropy;

}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){
    int area = rectArea(ul, lr);
    vector<int> dist = buildHist(ul, lr);
    return entropy(dist, area);
}
