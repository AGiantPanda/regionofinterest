#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "CImg.h"

/*struct POINT{
	int x;
	int y;
};*/

class MYPOINT{
public:
	MYPOINT() :x(0), y(0){}
	MYPOINT(int xx) :x(xx), y(xx){}
	MYPOINT(int xx, int yy) :x(xx), y(yy){}
	int x;//代表第几列
	int y;//代表第几行
};

//画线算法，在(x0, y0)和(x1, y1)之间画一条线并将线上的点填充成白色
void DrawLine(cimg_library::CImg<unsigned char>& mask, int x0, int y0, int x1, int y1){
	int dx = x1 - x0;
	int dy = y1 - y0;
	int ux = ((dx > 0) << 1) - 1;
	int uy = ((dy > 0) << 1) - 1;
	int x = x0, y = y0, eps;

	eps = 0; dx = abs(dx); dy = abs(dy);
	if (dx > dy){
		for (x = x0; x != x1; x += ux){
			mask.atXY(x,y) = 255;
			eps += dy;
			if ((eps << 1) >= dx){
				y += uy; eps -= dx;
			}
		}
	}
	else{
		for (y = y0; y != y1; y += uy){
			mask.atXY(x,y) = 255;
			eps += dx;
			if ((eps << 1) >= dy){
				x += ux; eps -= dy;
			}
		}
	}
}

/**
to judge whether a point is in a circle or not with acos(), return true/false
parameter:
@p		--the given point
@points	--the circle
note:
仅适用于凸多边形
*/
bool JudgePoint(MYPOINT p, std::vector<MYPOINT> Points)
{
	double angle = 0.0, dot;
	double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
	std::vector<MYPOINT>::iterator it1 = Points.begin();
	std::vector<MYPOINT>::iterator it2;
	for (it2 = (it1 + 1); it2 < Points.end(); ++it1, ++it2){
		x1 = (*it1).x - p.x;
		y1 = (*it1).y - p.y;
		x2 = (*it2).x - p.x;
		y2 = (*it2).y - p.y;
		dot = x1*x2 + y1*y2;
		dot = dot / (sqrt(x1*x1 + y1*y1)*sqrt(x2*x2 + y2*y2));
		angle += acos(dot) * 180 / 3.1415926;
	}
	it2 = Points.begin();
	x1 = (*it1).x - p.x;
	y1 = (*it1).y - p.y;
	x2 = (*it2).x - p.x;
	y2 = (*it2).y - p.y;
	dot = x1*x2 + y1*y2;
	dot = dot / (sqrt(x1*x1 + y1*y1)*sqrt(x2*x2 + y2*y2));
	angle += acos(dot) * 180 / 3.1415926;
	if (angle >= 360.0)
		return true;
	else
		return false;
}

//another
bool JudgePoint2(MYPOINT p, std::vector<MYPOINT> Points)
{
	if (Points.size() < 3)
		return false;
	int sum = 0;
	int count = Points.size();

	int x0, y0, x1, y1;
	int x;
	for (int i = 0; i < count; i++){
		if (i == count - 1){
			x0 = Points[i].x;
			y0 = Points[i].y;
			x1 = Points[0].x;
			y1 = Points[0].y;
		}
		else{
			x0 = Points[i].x;
			y0 = Points[i].y;
			x1 = Points[i + 1].x;
			y1 = Points[i + 1].y;
		}
		if ((p.y >= y0 && p.y < y1) || (p.y >= y1 && p.y < y0)){
			if (abs(y0 - y1)>0){
				x = x0 - (x0 - x1)*(y0 - p.y) / (y0 - y1);
				if (x < p.x)
					sum++;
			}
		}
	}
	if (sum % 2)
		return true;
	else 
		return false;
	return false;
}

//区域蒙版，为选定的区域做一个二值图像的蒙版，方便以后的计算
cimg_library::CImg<unsigned char> AreaMask(cimg_library::CImg<unsigned char> origin, std::vector<MYPOINT> Points){
	int i, j, x0 = origin.width(), y0 = origin.height(), x1 = 0, y1 = 0;
	cimg_library::CImg<unsigned char> mask(origin.width(), origin.height(), 1, 1);
	mask.fill(0);
	//draw the circle of chosen area
	for (i = 0, j = Points.size() - 1; i < Points.size(); j = i++){
		DrawLine(mask, Points[j].x, Points[j].y, Points[i].x, Points[i].y);
	}
	//draw the mask
	//find a point in the circle
	//缩小计算范围，找到能框住选定区域的方框，再从方框中部开始搜索在图像内部的点
	for (std::vector<MYPOINT>::iterator it = Points.begin(); it != Points.end(); it++){
		x0 = x0 < (*it).x ? x0 : (*it).x;
		y0 = y0 < (*it).y ? y0 : (*it).y;
		x1 = x1 > (*it).x ? x1 : (*it).x;
		y1 = y1 > (*it).y ? y1 : (*it).y;
	}
	MYPOINT p(x0, (y1+y0)/2);
	//while (!JudgePoint(p, Points) && p.y <= y1){
	//	p.x = (p.x + 1) % (x1 + 1) + x0;
	//	if (p.x == x0) p.y++;
	//}//cant use this function cause it takes too much computation
	int c = 0, tmp_x = Points[0].x, tmp_y = Points[0].y;
	while (p.y <= y1){
		if (mask.atXY(p.x, p.y) == 255) {
			if (c == 0){
				tmp_x = p.x + 1; tmp_y = p.y;
			}
			c++;
		}
		if (c == 2){
			if (mask.atXY(tmp_x, tmp_y) == 0 && tmp_y == p.y && tmp_x < p.x)
				break;
		}
		p.x = (p.x + 1) % (x1 + 1) + x0;
		if (p.x == x0){
			p.y++; c = 0;
		}
	}
	p.x = tmp_x; p.y = tmp_y;
	//std::cout << x0 << " " << y0 << std::endl;
	//std::cout << tmp_x << " " << tmp_y << std::endl;
	std::vector<MYPOINT> stack;
	stack.push_back(p);
	MYPOINT Cur = p, Las;
	while (stack.size()){
		mask.atXY(Cur.x, Cur.y) = 255;
		
		if (Cur.x + 1 < origin.width() && mask.atXY(Cur.x + 1, Cur.y) == 0){
			Cur.x += 1; stack.push_back(Cur); continue;
		}
		if (Cur.y + 1 < origin.height() && mask.atXY(Cur.x, Cur.y +1 ) == 0){
			Cur.y += 1; stack.push_back(Cur); continue;
		}
		if (Cur.x - 1 > -1 && mask.atXY(Cur.x - 1, Cur.y) == 0){
			Cur.x -= 1; stack.push_back(Cur); continue;
		}
		if (Cur.y - 1 > -1 && mask.atXY(Cur.x, Cur.y - 1) == 0){
			Cur.y -= 1; stack.push_back(Cur); continue;
		}
		stack.pop_back();
		if (stack.size()){
			Cur.x = stack.back().x;
			Cur.y = stack.back().y;
		}
	}
	return mask;
}

/**
otsu, return the computed threshold,the CImg version
parameter:
@origin		--original pic
@Points		--ROI, 给定连续的点，按顺序连城线以确定一个封闭的环
@debug		--debug option, 0 means no information outputed
*/
int otsu(cimg_library::CImg<unsigned char> origin, std::vector<MYPOINT> Points, cimg_library::CImg<unsigned char>& mask, int debug = 0){
	int thresholdValue = 1;
	int ihist[256] = { 0 };

	int i, j, k;
	int rows = origin.height(), cols = origin.width();
	int n, n1, n2, gmin = 255, gmax = 0;
	double m1, m2, sum, csum, fmax, sb;
	int x0 = origin.width(), y0 = origin.height(), x1 = 0, y1 = 0;
	//缩小计算范围
	for (std::vector<MYPOINT>::iterator it = Points.begin(); it != Points.end(); it++){
		x0 = x0 < (*it).x ? x0 : (*it).x;
		y0 = y0 < (*it).y ? y0 : (*it).y;
		x1 = x1 > (*it).x ? x1 : (*it).x;
		y1 = y1 > (*it).y ? y1 : (*it).y;
	}

	mask = AreaMask(origin, Points);
	

	for (i = y0; i <= y1; i++){
		for (j = x0; j < x1; j++){
			if (mask.atXY(j, i) == 255){
				ihist[origin.atXY(j,i)]++;
				if (origin.atXY(j,i) > gmax)gmax = origin.atXY(j,i);
				if (origin.atXY(j,i) < gmin)gmin = origin.atXY(j,i);
			}
		}
	}

	sum = csum = 0.0;
	n = 0;
	for (k = 0; k < 256; k++){
		sum += (double)k*(double)ihist[k];
		n += ihist[k];
	}

	if (!n){
		std::cout << "NOT NORMAL!!!" << std::endl << " thresholdValue = 160\n" << std::endl;
		return 160;
	}

	//do the otsu global thresholding method
	fmax = -1.0;
	n1 = 0;
	for (k = 0; k < 255; k++){
		n1 += ihist[k];
		if (!n1){ continue; }
		n2 = n - n1;
		if (n2 == 0){ break; }
		csum += (double)k*ihist[k];
		m1 = csum / n1;
		m2 = (sum - csum) / n2;
		sb = (double)n1*(double)n2*(m1 - m2)*(m1 - m2);
		if (sb>fmax){
			fmax = sb;
			thresholdValue = k;
		}
	}
	if (debug) std::cout << "#OTSU: thresholdValue = " << thresholdValue << std::endl
		<< "gmin = " << gmin << std::endl
		<< "gmax = " << gmax << std::endl;
	return thresholdValue;
}

/**
find all the runs in the img
parameter:
@origin				--the original CImg pic
@NumberofRuns		--count the number of runs
@stRun				--record the start of a run
@enRun				--record the end of a run
@rowRun				--record which row the run is in
@masRun				--record how much mass a run has
*/
void fillRunVector(cimg_library::CImg<unsigned char> origin, int &NumberofRuns, std::vector<int> &stRun, std::vector<int> &enRun, std::vector<int> &rowRun, std::vector<int> &masRun)
{
	for (int i = 0; i < origin.height(); i++){
		if (origin.atXY(0, i) == 255){
			NumberofRuns++;
			stRun.push_back(0);
			rowRun.push_back(i);
		}
		for (int j = 1; j < origin.width(); j++){
			if (origin.atXY(j - 1, i) == 0 && origin.atXY(j, i) == 255){
				NumberofRuns++;
				stRun.push_back(j);
				rowRun.push_back(i);
			}
			else if (origin.atXY(j-1, i) == 255 && origin.atXY(j, i) == 0) {
				enRun.push_back(j - 1);
				masRun.push_back(enRun[NumberofRuns-1] - stRun[NumberofRuns-1] + 1);
			}
		}
		if (origin.atXY(origin.width() - 1,i)){
			enRun.push_back(origin.width() - 1);
			masRun.push_back(enRun[NumberofRuns-1] - stRun[NumberofRuns-1] + 1);
		}
	}
}

/**
完成团的标记以及等价对列表的生成
对所有团进行标记，并将上下两行有重叠的团标记到等价序列中
@stRun, enRun, rowRun, NumberofRuns		--runs
@runLabels	equivalences				--label & equivalences
@offset									--8连通域或是4连通域
*/
void firstPass(std::vector<int>& stRun, std::vector<int>& enRun, std::vector<int>& rowRun, int NumberofRuns,
	std::vector<int>& runLabels, std::vector<std::pair<int, int>>& equivalences, int offset)
{
	runLabels.assign(NumberofRuns, 0);
	int idxLabel = 1;
	int curRowIdx = 0;
	int firstRunOnCur = 0;
	int firstRunOnPre = 0;
	int lastRunOnPre = -1;
	for (int i = 0; i < NumberofRuns; i++){
		if (rowRun[i] != curRowIdx){
			curRowIdx = rowRun[i];
			firstRunOnPre = firstRunOnCur;
			lastRunOnPre = i - 1;
			firstRunOnCur = i;
		}
		for (int j = firstRunOnPre; j <= lastRunOnPre; j++){
			if (stRun[i] <= enRun[j] + offset && enRun[i] >= stRun[j] - offset && rowRun[i] - rowRun[j] == 1){
				if (runLabels[i] == 0)//没有被标记
					runLabels[i] = runLabels[j];
				else if (runLabels[i] != runLabels[j])//已经被标记
					equivalences.push_back(std::make_pair(runLabels[i], runLabels[j]));//加入等价对
			}
		}
		if (runLabels[i] == 0)//没有与上一行重合的run
			runLabels[i] = idxLabel++;
	}
}

/**
等价对的处理，将等价对转化为若干个等价序列，
这个过程中将每个团都看成一个图的节点，等价对表示两个节点之间有通路，由此形成无向图
采用深度优先的遍历算法，找到所有的序列
每个等价表用一个vector<int>来保存，等价对列表保存在map<pair<int, int>中
*/
void replaceSameLabel(std::vector<int>& runLabels, std::vector<std::pair<int, int>>& equivalence)
{
	int maxLabel = *std::max_element(runLabels.begin(), runLabels.end());
	std::vector<std::vector<bool>> eqTab(maxLabel, std::vector<bool>(maxLabel, false));//maxLabel*maxLabel的bool型二维数组，值全为false
	std::vector<std::pair<int, int>>::iterator vecPairIt = equivalence.begin();
	while (vecPairIt != equivalence.end()){
		eqTab[vecPairIt->first - 1][vecPairIt->second - 1] = true;
		eqTab[vecPairIt->second - 1][vecPairIt->first - 1] = true;
		vecPairIt++;
	}
	std::vector<int> labelFlag(maxLabel, 0);
	std::vector<std::vector<int>> equaList;
	std::vector<int> tempList;
	//std::cout << maxLabel << std::endl;
	for (int i = 1; i <= maxLabel; i++){
		if (labelFlag[i - 1]){
			continue;
		}
		labelFlag[i - 1] = equaList.size() + 1;
		tempList.push_back(i);
		for (std::vector<int>::size_type j = 0; j < tempList.size(); j++){
			for (std::vector<bool>::size_type k = 0; k != eqTab[tempList[j] - 1].size(); k++){
				if (eqTab[tempList[j] - 1][k] && !labelFlag[k]){
					tempList.push_back(k + 1);
					labelFlag[k] = equaList.size() + 1;
				}
			}
		}
		equaList.push_back(tempList);
		tempList.clear();
	}
	//std::cout << equaList.size() << std::endl;
	for (std::vector<int>::size_type i = 0; i != runLabels.size(); i++){
		runLabels[i] = labelFlag[runLabels[i] - 1];
	}
}

/*
输入原二值图像，返回一个MYPOINT类型的容器，容器内标记最大连通域的外轮廓
这是直接对已经找到最大连通域的二值图像进行外轮廓标记
*/
std::vector<MYPOINT> FindContour(cimg_library::CImg<unsigned char> origin)
{
	int i = 0, j = 0;
	int cols = origin.width(), rows = origin.height();
	int imgsize = (cols + 2)*(rows + 2);
	unsigned char *ori = new unsigned char[imgsize];
	memset(ori, 0, imgsize*sizeof(unsigned char));
	for (i = 0; i < rows; i++){
		for (j = 0; j < cols; j++)
			ori[(i + 1)*(cols + 2) + (j + 1)] = origin[j,i];
	}
	int *oriLabel = new int[imgsize];
	memset(oriLabel, 0, imgsize*sizeof(int));
	std::vector < MYPOINT > contour;

	//find the start point
	for (i = 0; i < rows + 2; i++){
		for (j = 0; j < cols + 2 && ori[i*(cols + 2) + j] != 255; j++){}
		if (ori[i*(cols + 2) + j] == 255)
			break;
	}

	if (i == rows + 2 && j == cols + 2 && ori[imgsize] == 0) {
		std::cout << "there is no contour." << std::endl;
		return contour;
	}
	MYPOINT start(i - 1, j - 1);
	contour.push_back(start);
	oriLabel[i*(cols + 2) + j] = 1;

	//设置循环计算的参数
	MYPOINT CurP;
	MYPOINT LastP = start;
	int x[8] = { i, i + 1, i + 1, i + 1, i, i - 1, i - 1, i - 1 };
	int y[8] = { j + 1, j + 1, j, j - 1, j - 1, j - 1, j, j + 1 };
	int c = 7, n = 0;
	for (n = 0; n < 8; n++){
		if (ori[x[c] * (cols + 2) + y[c]] == 255) {
			oriLabel[x[c] * (cols + 2) + y[c]] = 1;
			break;
		}
		else
			oriLabel[x[c] * (cols + 2) + y[c]] = -1;
		c = (c + 1) % 8;
	}
	if (n == 7 && ori[x[c] * (cols + 2) + y[c]] == 0) {
		return contour;
	}
	else{
		CurP.x = x[c] - 1;
		CurP.y = y[c] - 1;
		contour.push_back(CurP);
		//重新计算x[]，y[]，c， n
		x[0] = CurP.x + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1 + 1, x[3] = CurP.x + 1 + 1, x[4] = CurP.x + 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1 - 1, x[7] = CurP.x + 1 - 1;
		y[0] = CurP.y + 1 + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1, y[3] = CurP.y + 1 - 1, y[4] = CurP.y + 1 - 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1, y[7] = CurP.y + 1 + 1;
		n = 0, c = (c + 6) % 8;
	}
	//std::cout << start << std::endl;
	int stop = 1;
	while (stop != 0){
		for (n = 0; n < 8; n++){
			if (ori[x[c] * (cols + 2) + y[c]] == 255) {
				oriLabel[x[c] * (cols + 2) + y[c]] = 1;
				break;
			}
			else
				oriLabel[x[c] * (cols + 2) + y[c]] = -1;
			c = (c + 1) % 8;
		}
		LastP = CurP;
		CurP.x = x[c] - 1;
		CurP.y = y[c] - 1;
		contour.push_back(CurP);
		//重新计算x[]，y[]，c， n
		if (LastP.x == start.x && LastP.y == start.y){
			//for (int cc = 0; cc < 7; cc++)
			//	printf("%3d, %3d, %3d, %3d, %3d, %3d, %3d\n", oriLabel[cc * 7 + 0], oriLabel[cc * 7 + 1], oriLabel[cc * 7 + 2], oriLabel[cc * 7 + 3], oriLabel[cc * 7 + 4], oriLabel[cc * 7 + 5], oriLabel[cc * 7 + 6], oriLabel[cc * 7 + 7]);
			for (n = 0; n < 8 && (oriLabel[x[n] * (cols + 2) + y[n]] != 0 || ori[x[n] * (cols + 2) + y[n]] == 255); n++);
			if (n == 8 && (oriLabel[x[7] * (cols + 2) + y[7]] != 0 || ori[x[7] * (cols + 2) + y[7]] == 255))
				stop = 0;
			contour.pop_back();
		}
		x[0] = CurP.x + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1 + 1, x[3] = CurP.x + 1 + 1, x[4] = CurP.x + 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1 - 1, x[7] = CurP.x + 1 - 1;
		y[0] = CurP.y + 1 + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1, y[3] = CurP.y + 1 - 1, y[4] = CurP.y + 1 - 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1, y[7] = CurP.y + 1 + 1;
		c = (c + 6) % 8;
	}

	//unsigned char *op = new unsigned char[rows*cols];
	//for (int i = 0; i < rows; i++){
	//	for (int j = 0; j < cols; j++){
	//		if (oriLabel[(i + 1)*(cols + 2) + (j + 1)] == 1)
	//			op[i*cols + j] = 255;
	//		else
	//			op[i*cols + j] = 0;
	//	}
	//}
	//cimg_library::CImg<unsigned char> ttt(op, rows, cols);
	//cimg_library::CImgDisplay ddd(ttt, "ddd");
	return contour;
}

std::vector<MYPOINT> FindContour(const unsigned char *origin, int rows, int cols)
{
	int i = 0, j = 0;
	int imgsize = (cols + 2)*(rows + 2);
	unsigned char *ori = new unsigned char[imgsize];
	memset(ori, 0, imgsize*sizeof(unsigned char));
	for (i = 0; i < rows; i++){
		for (j = 0; j < cols; j++)
			ori[(i + 1)*(cols + 2) + (j + 1)] = origin[i*cols + j];
	}
	int *oriLabel = new int[imgsize];
	memset(oriLabel, 0, imgsize*sizeof(int));
	std::vector < MYPOINT > contour;

	//find the start MYPOINT
	for (i = 0; i < rows + 2; i++){
		for (j = 0; j < cols + 2 && ori[i*(cols + 2) + j] != 255; j++){}
		if (ori[i*(cols + 2) + j] == 255)
			break;
	}

	if (i == rows + 2 && j == cols + 2 && ori[imgsize] == 0) {
		std::cout << "there is no contour." << std::endl;
		return contour;
	}
	MYPOINT start(i - 1, j - 1);
	contour.push_back(start);
	oriLabel[i*(cols + 2) + j] = 1;

	//设置循环计算的参数
	MYPOINT CurP;
	MYPOINT LastP = start;
	int x[8] = { i, i + 1, i + 1, i + 1, i, i - 1, i - 1, i - 1 };
	int y[8] = { j + 1, j + 1, j, j - 1, j - 1, j - 1, j, j + 1 };
	int c = 7, n = 0;
	for (n = 0; n < 8; n++){
		if (ori[x[c] * (cols + 2) + y[c]] == 255) {
			oriLabel[x[c] * (cols + 2) + y[c]] = 1;
			break;
		}
		else
			oriLabel[x[c] * (cols + 2) + y[c]] = -1;
		c = (c + 1) % 8;
	}
	if (n == 7 && ori[x[c] * (cols + 2) + y[c]] == 0) {
		return contour;
	}
	else{
		CurP.x = x[c] - 1;
		CurP.y = y[c] - 1;
		contour.push_back(MYPOINT(CurP.y,CurP.x));
		//重新计算x[]，y[]，c， n
		x[0] = CurP.x + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1 + 1, x[3] = CurP.x + 1 + 1, x[4] = CurP.x + 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1 - 1, x[7] = CurP.x + 1 - 1;
		y[0] = CurP.y + 1 + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1, y[3] = CurP.y + 1 - 1, y[4] = CurP.y + 1 - 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1, y[7] = CurP.y + 1 + 1;
		n = 0, c = (c + 6) % 8;
	}
	//std::cout << start << std::endl;
	int stop = 1;
	while (stop != 0){
		for (n = 0; n < 8; n++){
			if (ori[x[c] * (cols + 2) + y[c]] == 255) {
				oriLabel[x[c] * (cols + 2) + y[c]] = 1;
				break;
			}
			else
				oriLabel[x[c] * (cols + 2) + y[c]] = -1;
			c = (c + 1) % 8;
		}
		LastP = CurP;
		CurP.x = x[c] - 1;
		CurP.y = y[c] - 1;
		contour.push_back(MYPOINT(CurP.y, CurP.x));
		//重新计算x[]，y[]，c， n
		if (LastP.x == start.x && LastP.y == start.y){
			//for (int cc = 0; cc < 7; cc++)
			//	printf("%3d, %3d, %3d, %3d, %3d, %3d, %3d\n", oriLabel[cc * 7 + 0], oriLabel[cc * 7 + 1], oriLabel[cc * 7 + 2], oriLabel[cc * 7 + 3], oriLabel[cc * 7 + 4], oriLabel[cc * 7 + 5], oriLabel[cc * 7 + 6], oriLabel[cc * 7 + 7]);
			for (n = 0; n < 8 && (oriLabel[x[n] * (cols + 2) + y[n]] != 0 || ori[x[n] * (cols + 2) + y[n]] == 255); n++);
			if (n == 8 && (oriLabel[x[7] * (cols + 2) + y[7]] != 0 || ori[x[7] * (cols + 2) + y[7]] == 255))
				stop = 0;
			contour.pop_back();
		}
		x[0] = CurP.x + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1 + 1, x[3] = CurP.x + 1 + 1, x[4] = CurP.x + 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1 - 1, x[7] = CurP.x + 1 - 1;
		y[0] = CurP.y + 1 + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1, y[3] = CurP.y + 1 - 1, y[4] = CurP.y + 1 - 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1, y[7] = CurP.y + 1 + 1;
		c = (c + 6) % 8;
	}

	return contour;
}

//修改一下上面的代码让其实现找所有团的外轮廓

//找最大的团并返回外轮廓
std::vector<MYPOINT> FindBiggestContour(cimg_library::CImg<unsigned char> origin, std::vector<MYPOINT> Points)
{
	std::vector<MYPOINT> contour;
	cimg_library::CImg<unsigned char> mask;
	//compute the binary img
	int threshold = otsu(origin, Points, mask);
	unsigned char *bi_img = new unsigned char[origin.width()*origin.height()];
	for (int i = 0; i < origin.height(); i++){
		for (int j = 0; j < origin.width(); j++){
			if (origin.atXY(j, i) >= threshold && mask.atXY(j, i) == 255)
				bi_img[i*origin.width() + j] = 255;
			else
				bi_img[i*origin.width() + j] = 0;
		}
	}
	cimg_library::CImg<unsigned char> bi(bi_img, origin.width(), origin.height());
	//done

	std::vector<int> stRun, enRun, rowRun, masRun;
	int NumberofRuns = 0, offset = 1, maxRun = 1;
	fillRunVector(bi, NumberofRuns, stRun, enRun, rowRun, masRun);
	std::vector<int> runLabels;
	std::vector<std::pair<int, int>> equivalences;
	firstPass(stRun, enRun, rowRun, NumberofRuns, runLabels, equivalences, offset);
	if (!NumberofRuns){
		std::cout << "NOTHING FOUND!!!" << std::endl 
			<<"size of contour is: "<< contour.size() <<  std::endl;
		return contour;
	}
	replaceSameLabel(runLabels, equivalences);

	int maxLabel = *max_element(runLabels.begin(), runLabels.end());
	int *MassofRuns = new int[maxLabel + 1];
	memset(MassofRuns, 0, (maxLabel + 1)*sizeof(int));

	memset(bi_img, 0, (bi.height()*bi.width())*sizeof(unsigned char));
	for (int c = 0; c < NumberofRuns; c++){
		MassofRuns[runLabels[c]] += masRun[c];
		if (MassofRuns[runLabels[c]]>maxRun){
			maxRun = MassofRuns[runLabels[c]];
			maxLabel = runLabels[c];
		}
		int i = rowRun[c]; 
		for (int j = stRun[c]; j <= enRun[c]; j++){
			bi_img[i*bi.width() + j] = runLabels[c];
		}
	}
	for (int i = 0; i < bi.height(); i++)
	for (int j = 0; j < bi.width(); j++) {
		if (bi_img[i*bi.width()+j] == maxLabel)
			bi_img[i*bi.width() + j] = 255;
		else
			bi_img[i*bi.width() + j] = 0;
	}

	return FindContour(bi_img, bi.width(), bi.height());
}

