
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
	:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
	{}

toqutree::~toqutree(){
	clear(root);
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
	}
	return *this;
}

toqutree::toqutree(PNG & imIn, int k){ 

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */

/* your code here */
	int centerX = imIn.width()/2;
	int centerY = imIn.height()/2;
	int ulX = centerX - pow(2, k)/2;
	int ulY = centerY - pow(2, k)/2;
	PNG* sub = new PNG(pow(2, k),pow(2, k));
	for (int i = 0; i < pow(2, k); i++){
		for (int j = 0; j < pow(2, k); j++){
			int x = i + ulX;
			int y = j + ulY;
			HSLAPixel* pixel = imIn.getPixel(x, y);
			HSLAPixel* newPixel = sub->getPixel(i, j);
			*newPixel = *pixel;
		}
	}
	//printf("%d\n", 55);
	root = buildTree(sub, k);
}

int toqutree::size() {
/* your code here */
	return size(root);
}

int toqutree::size(const Node* croot){
	if(croot != NULL){
		return size(croot->SE) + size(croot->SW) + size(croot->NE) + size(croot->NW) + 1;
	}else
		return 0;
}
int debug = 0;

toqutree::Node * toqutree::buildTree(PNG * im, int k) {

/* your code here */

// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used .
// similarly, at each level of the tree you will want to 
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.
	stats* s = new stats(*im);

	pair<int, int> ul(0, 0);
	pair<int, int> lr((int)im->width() - 1, (int)im->height() - 1);

	//PNG* empty = new PNG((int)im->width(), (int)im->height());
	//if(*im == *empty)
		//printf("%d\n", 9999999);

	PNG* SE;
	PNG* SW;
	PNG* NE;
	PNG* NW;
	
	HSLAPixel avg = s->getAvg(ul, lr);
	//printf("%f\n", avg.h);

	Node* node;

	if (k == 0){
		node = new Node(ul, k, *im->getPixel(0,0));

		delete s;
		delete im;
		return node;
	}
	if (k == 1){
		node = new Node(lr, k, avg);

		SE = new PNG(1,1);
		SW = new PNG(1,1);
		NE = new PNG(1,1);
		NW = new PNG(1,1);

		*SE->getPixel(0,0) = *im->getPixel(1,1);
		*SW->getPixel(0,0) = *im->getPixel(0,1);
		*NE->getPixel(0,0) = *im->getPixel(1,0);
		*NW->getPixel(0,0) = *im->getPixel(0,0);

		node->SE = buildTree(SE, k-1);
		node->SW = buildTree(SW, k-1);
		node->NE = buildTree(NE, k-1);
		node->NW = buildTree(NW, k-1);

		delete s;
		delete im;
		return node;
	}
	SE = new PNG(pow(2, k-1), pow(2, k-1));
	SW = new PNG(pow(2, k-1), pow(2, k-1));
	NE = new PNG(pow(2, k-1), pow(2, k-1));
	NW = new PNG(pow(2, k-1), pow(2, k-1));

	pair<int, int> newUl(ul.first + pow(2, k-1)/2, ul.second + pow(2, k-1)/2);
	pair<int, int> newLr(lr.first - pow(2, k-1)/2, lr.second - pow(2, k-1)/2);

	int area = (pow(2,k-1) * pow(2,k-1));
		
	pair<int, int> ctr(0, 0);

	double min_entropy = 10000000000;

	for (int i = newUl.first; i < newLr.first+1; i++){
		for (int j = newUl.second; j < newLr.second+1; j++){
			double newEntropy;
			pair<int, int> newCtr(i, j);
			vector<int> seDist;
			vector<int> swDist;
			vector<int> neDist;
			vector<int> nwDist;
			if((i + pow(2, k-1) - 1) < lr.first && (j + pow(2, k-1) - 1) < lr.second){
				pair<int, int> seLr((i + pow(2, k-1) - 1),  (j + pow(2, k-1) - 1));

				seDist = s->buildHist(newCtr, seLr);
				swDist.resize(36);
				neDist.resize(36);
				nwDist.resize(36);
				for (int h = 0; h < 36; h++){
					swDist[h] = s->hist[lr.first][j + pow(2, k-1) - 1][h] - s->hist[lr.first][j - 1][h] - seDist[h];
				}
				for (int h = 0; h < 36; h++){
					neDist[h] = s->hist[i + pow(2, k-1) - 1][lr.second][h] - s->hist[i - 1][lr.second][h] - seDist[h];
				}
				for (int h = 0; h < 36; h++){
					nwDist[h] = s->hist[lr.first][lr.second][h] - neDist[h] - swDist[h] - seDist[h];
				}		
			}else if((i + pow(2, k-1) - 1) < lr.first && (j + pow(2, k-1) - 1) >= lr.second){
				pair<int, int> neUl(i, j - pow(2, k - 1));
				pair<int, int> neLr((i + pow(2, k - 1) - 1), j - 1);

				neDist = s->buildHist(neUl, neLr);
				swDist.resize(36);
				seDist.resize(36);
				nwDist.resize(36);
				for (int h = 0; h < 36; h++){
					if((j + pow(2, k-1) - 1) == lr.second)
						nwDist[h] = s->hist[lr.first][j - 1][h] - neDist[h];
					else
						nwDist[h] = s->hist[lr.first][j - 1][h] - s->hist[lr.first][neUl.second - 1][h] - neDist[h];
					
				}
				for (int h = 0; h < 36; h++){
					seDist[h] = s->hist[i + pow(2, k-1) - 1][lr.second][h] - s->hist[i - 1][lr.second][h] - neDist[h];
				}
				for (int h = 0; h < 36; h++){
					swDist[h] = s->hist[lr.first][lr.second][h] - neDist[h] - nwDist[h] - seDist[h];
				}
			
			}else if((i + pow(2, k-1) - 1) >= lr.first && (j + pow(2, k-1) - 1) < lr.second){
				pair<int, int> swUl((i - pow(2, k-1)), j);
				pair<int, int> swLr(i - 1, (j + pow(2, k-1) - 1));

				swDist = s->buildHist(swUl, swLr);
				seDist.resize(36);
				neDist.resize(36);
				nwDist.resize(36);
				for (int h = 0; h < 36; h++){
					seDist[h] = s->hist[lr.first][j + pow(2, k-1) - 1][h] - s->hist[lr.first][j - 1][h] - swDist[h];
				}
				for (int h = 0; h < 36; h++){
					if((i + pow(2, k-1) - 1) == lr.first)
						nwDist[h] = s->hist[i - 1][lr.second][h] - swDist[h];
					else
						nwDist[h] = s->hist[i - 1][lr.second][h] - s->hist[swUl.first - 1][lr.second][h] - swDist[h];
				
				}
				for (int h = 0; h < 36; h++){
					neDist[h] = s->hist[lr.first][lr.second][h] - swDist[h] - nwDist[h] - seDist[h];
				}
			}else{
				pair<int, int> nwUl((i - pow(2, k-1)), (j - pow(2, k-1)));
				pair<int, int> nwLr(i - 1, j - 1);

				nwDist = s->buildHist(nwUl, nwLr);
				swDist.resize(36);
				neDist.resize(36);
				seDist.resize(36);
				for (int h = 0; h < 36; h++){
					if((j + pow(2, k-1) - 1) == lr.second)
						neDist[h] = s->hist[lr.first][j - 1][h] - nwDist[h];
					else
						neDist[h] = s->hist[lr.first][j - 1][h] - s->hist[lr.first][nwUl.second - 1][h] - nwDist[h];
				}
				for (int h = 0; h < 36; h++){
					if((i + pow(2, k-1) - 1) == lr.first)
						swDist[h] = s->hist[i - 1][lr.second][h] - nwDist[h];
					else 
						swDist[h] = s->hist[i - 1][lr.second][h] - s->hist[nwUl.first - 1][lr.second][h] - nwDist[h];
				}
				for (int h = 0; h < 36; h++){
					seDist[h] = s->hist[lr.first][lr.second][h] - swDist[h] - nwDist[h] - neDist[h];
				}
			}
			newEntropy = s->entropy(seDist, area) + s->entropy(swDist, area) + s->entropy(neDist, area) + s->entropy(nwDist, area);
			
			if (newEntropy < min_entropy){
				min_entropy = newEntropy;
				ctr.first = newCtr.first;
				ctr.second = newCtr.second;
			}
		} 
	}

	node = new Node(ctr, k, avg);

	int newL = pow(2,k-1);
    int SE_X = ctr.first;
    int SE_Y = ctr.second;
	int sl = pow(2,k);
    //debug++;
    //if(debug>8 ) return NULL;
    //printf("Debug:%d entropy:%f Origlength:%d x:%d  y:%d\n",debug,min_entropy,sl,SE_X, SE_Y);

    // set up SE 
    //SE = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+i)%sl;
            int y = (SE_Y+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = SE->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }

    // set up NE
    //NE = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+i)%sl;
            int y = (SE_Y+newL+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = NE->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }

    // set up SW
    //SW = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+newL+i)%sl;
            int y = (SE_Y+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = SW->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }

    // set up NW
    //NW = new PNG(newL,newL);
    for(int i = 0; i < newL; i++){
        for (int j = 0; j < newL; j++){
            int x = (SE_X+newL+i)%sl;
            int y = (SE_Y+newL+j)%sl;
            HSLAPixel* p = im->getPixel(x,y);
            HSLAPixel* c = NW->getPixel(i,j);
            c->h = p->h;
            c->s = p->s;
            c->l = p->l;
            c->a = p->a;
        }
    }

    // printf("%f\n", min_entropy);

    // croot = new Node(split_point,k,im_stats->getAvg(make_pair(0,0), make_pair((int)im->width() - 1, (int)im->height() - 1)));

    node->SE = buildTree(SE, k-1);
    node->NW = buildTree(NW, k-1);
    node->NE = buildTree(NE, k-1);
    node->SW = buildTree(SW, k-1);

	delete s;
    delete im;
    return node;
	// node->SE = buildTree(SE, k-1);
	// node->SW = buildTree(SW, k-1);
	// node->NE = buildTree(NE, k-1);
	// node->NW = buildTree(NW, k-1);
	
	// delete im;
	// return node;
}

PNG toqutree::render(){

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the 
// quadtree, instead.

/* your code here */
	PNG png; png.resize(pow(2, root->dimension), pow(2, root->dimension));
	for (int i = 0; i < (int)png.width(); i++){
		for (int j = 0; j < (int)png.height(); j++){
			pair<int, int> cord(i, j);
			Node* node = find(root, cord);
			*png.getPixel(i, j) = node->avg;
		}
	}
	return png;
}

toqutree::Node * toqutree::find(Node* root, const pair<int, int> &key){
    // todo
    pair<int, int> ctr = root->center;
	int k = root->dimension;
    int origRend = pow(2,k) - 1;
    int subLength = pow(2, k - 1);
	if (root->SE == NULL || k == 0) return root;

	else if(((ctr.first + subLength) > origRend) && ((ctr.second + subLength) <=  origRend)){
		if((((key.first >= ctr.first) && (key.first <= pow(2, k) - 1)) || ((key.first >= 0) && (key.first < ctr.first-subLength))) && 
            (key.second >= ctr.second) && (key.second < ctr.second+subLength)){
                int newX = ctr.first;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->SE, newKey);
            }
		else if((((key.second >= ctr.second + subLength) && (key.second <= pow(2, k) - 1)) || ((key.second>=0) && (key.second < ctr.second))) &&
            (key.first >= ctr.first - subLength) && (key.first < ctr.first)){
                int newX = ctr.first - subLength;
                int newY = ctr.second + subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->NW, newKey);
            }
		else if((key.second >= ctr.second) && (key.second < (ctr.second + subLength)) && (key.first < ctr.first) && (key.first >= (ctr.first - subLength))){
            int newX = ctr.first - subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->SW, newKey);
        }
		else{
            int newX = ctr.first;
            int newY = ctr.second +subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->NE, newKey);
        }
			
	}
	else if(((ctr.first + subLength) <= origRend) && ((ctr.second + subLength) <= origRend)){
		if((key.first < (ctr.first + subLength)) && (key.second < (ctr.second + subLength)) && (key.first >= ctr.first) && (key.second >= ctr.second)){
            pair<int,int> newKey = make_pair(key.first-ctr.first, key.second-ctr.second);
            return find(root->SE, newKey);
        }
		else if((key.first >= ctr.first) && (key.first < (ctr.first + subLength)) && 
        (((key.second < ctr.second) && (key.second >=0)) || ((key.second< pow(2, k)) && (key.second >= (ctr.second + subLength))))){
			int newX = ctr.first;
            int newY = ctr.second+subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->NE, newKey);
        }
		else if((key.second >= ctr.second) && (key.second < (ctr.second + subLength)) && ((key.first < ctr.first) || (key.first >= (ctr.first + subLength)))){
			int newX = ctr.first + subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->SW, newKey);
        }
		else{
            int newX = ctr.first +subLength;
            int newY = ctr.second + subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->NW, newKey);
        }
	}
	else if(((ctr.first + subLength) <= origRend) && ((ctr.second + subLength) > pow(2, k) - 1)){
        if((key.first < (ctr.first + subLength)) && (key.first >= ctr.first) && 
            ((key.second >= ctr.second && key.second <= pow(2, k) - 1) || (key.second>=0 && key.second < (ctr.second - subLength)))){
                int newX = ctr.first;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->SE, newKey);
            }
        else if((key.first >= ctr.first) && (key.first < (ctr.first + subLength)) && (key.second < ctr.second) && (key.second >= (ctr.second - subLength))){
                int newX = ctr.first;
                int newY = ctr.second -subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->NE, newKey);
        }
        else if((key.second < ctr.second) && (key.second >= (ctr.second - subLength)) && 
            (((key.first >= ctr.first+subLength) && (key.first < pow(2, k))) || ((key.first>=0) && (key.first < ctr.first)))){
                int newX = ctr.first + subLength;
                int newY = ctr.second -subLength;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->NW, newKey);
            }
		else{
                int newX = ctr.first + subLength;
                int newY = ctr.second;
                pair<int,int> newKey = findNewKey(key, k, newX, newY);
                return find(root->SW, newKey);
        }
	}
    else {
		if((key.first >= (ctr.first - subLength)) && (key.second >= (ctr.second - subLength)) && (key.first < ctr.first) && (key.second < ctr.second)){
            int newX = ctr.first - subLength;
            int newY = ctr.second - subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->NW, newKey);
        }
		else if((((key.first >= ctr.first) && (key.first < pow(2, k))) || ((key.first >= 0) && (key.first < ctr.first-subLength))) && 
        (key.second >= ctr.second - subLength) && (key.second < ctr.second)){
            int newX = ctr.first;
            int newY = ctr.second -subLength;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->NE, newKey);
        }
		else if((key.first >= (ctr.first - subLength)) && (key.first < ctr.first)){
            int newX = ctr.first-subLength;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->SW, newKey);
        }
		else{
            int newX = ctr.first;
            int newY = ctr.second;
            pair<int,int> newKey = findNewKey(key, k, newX, newY);
            return find(root->SE, newKey);
        }
	}
}

pair<int,int> toqutree::findNewKey(pair<int,int> key, int k, int x, int y){
    int newF;
    int newS;
    
    if(key.first >= x){
        newF = key.first - x;
    }else{
        newF = key.first + (pow(2,k) - x);
    }

    if(key.second >= y){
        newS = key.second - y;
    }else{
        newS = key.second + (pow(2,k) - y);
    }

    return make_pair(newF, newS);
}

/* oops, i left the implementation of this one in the file! */
void toqutree::prune(double tol){

	prune(root,tol);

}

void toqutree::prune(Node* croot, double tol){
	if(croot != NULL){
		if(withTol(croot, tol, croot->avg)){
			clear(croot->SE);
			clear(croot->SW);
			clear(croot->NE);
			clear(croot->NW);
		}else{
		// int seH = hPrune(croot->SE, tol);
		// int swH = hPrune(croot->SW, tol);
		// int nwH = hPrune(croot->NW, tol);
		// int neH = hPrune(croot->NE, tol);

		// printf("%d\n", seH);
		// printf("%d\n", swH);
		// printf("%d\n", nwH);
		// printf("%d\n", neH);
		
		// if(seH > swH && seH > nwH && seH > neH)
		// 	prune(croot->SE, tol);
		// else if(swH > seH && swH > nwH && swH > neH)
		// 	prune(croot->SW, tol);
		// else if(neH > swH && neH > nwH && neH > seH)
		// 	prune(croot->NE, tol);
		// else
		// 	prune(croot->NW, tol);
			prune(croot->SE, tol);
			prune(croot->NE, tol);
			prune(croot->SW, tol);
			prune(croot->NW, tol);
		}

	}
}

int toqutree::hPrune(Node* croot, double tol){
	if (withTol(croot, tol, croot->avg)){
		return croot->dimension;
	}else{
		int max[] = {hPrune(croot->SE, tol), hPrune(croot->SW, tol), hPrune(croot->NE, tol), hPrune(croot->NW, tol)};
		return *max_element(max, max+4);
	}
}

bool toqutree::withTol(Node* croot, double tol, HSLAPixel avgPixel){
	if(croot->SE == NULL && croot->SW == NULL && croot->NE == NULL && croot->NW == NULL)
		return (avgPixel.dist(croot->avg) <= tol);
	else
		return withTol(croot->SE, tol, avgPixel) && withTol(croot->SW, tol, avgPixel) && withTol(croot->NE, tol, avgPixel) && withTol(croot->NW, tol, avgPixel);
	
	
}

/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr){
/* your code here */
	if (curr != NULL){
		clear(curr->SE);
		clear(curr->SW);
		clear(curr->NE);
		clear(curr->NW);
		delete curr;
		curr = NULL;
	}
}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {

/* your code here */
	Node* root = NULL;
	if(other != NULL){
		root = new Node(other->center, other->dimension, other->avg);
		root->SE = copy(other->SE);
		root->SW = copy(other->SW);
		root->NE = copy(other->NE);
		root->NW = copy(other->NW);
	}
	return root;
}


