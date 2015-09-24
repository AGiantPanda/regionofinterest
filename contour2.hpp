#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <list>
#include "CImg.h"

class Point{
public:
	Point() :x(0), y(0){}
	Point(int xx) :x(xx), y(xx){}
	Point(int xx, int yy) :x(xx), y(yy){}
	int x;//代表第几列
	int y;//代表第几行
};

template<typename T>
bool is_inPolygon(T x, T y, std::vector<Point> points)
{
	if (points.size() < 3){
		std::cout << "that is a line" << std::endl;
		return false;
	}

	int sum = 0;
	int count = points.size();

	float x0, y0, x1, y1;
	float tmp_x;
	for (int i = 0; i < count; i++){
		if (i == count - 1){
			x0 = float(points[i].x); y0 = float(points[i].y);
			x1 = float(points[0].x); y1 = float(points[0].y);
		}
		else{
			x0 = float(points[i].x); y0 = float(points[i].y);
			x1 = float(points[i + 1].x); y1 = float(points[i + 1].y);
		}
		if ((y >= y0 && y < y1) || (y >= y1 && y < y0)){
			if (abs(y0 - y1)>0){
				tmp_x = x0 - (x0 - x1)*(y0 - y) / (y0 - y1);
				if (tmp_x < x)
					sum++;
			}
		}
	}
	if (sum % 2)
		return true;
	else
		return false;
	return false;
}//for a point

template<typename T>
bool is_inPolygon(T x0, T y0, T x1, T y1, std::vector<Point> points)
{

	return true;
}//for a line

//前端传过来的是点的话，将这些点连线形成一个封闭图形
template<typename T>
void draw_line_loop(cimg_library::CImg<T>& mask, std::vector<Point> points)
{
	const unsigned char white[] = { 255, 255, 255 };

	//itrator the points and draw the lines
	for (std::vector<int>::size_type i = 0, j = points.size() - 1; i < points.size(); j = i++){
		mask.draw_line(points[i].x, points[i].y, points[j].x, points[j].y, white);
	}
}//pass

//4连通域的区域蒙版函数，根据提供的点返回一个同类型的蒙版CImg类，蒙版为1通道二值图像。
//速度慢，应该是在容器的操作上浪费了比较多的时间
template<typename T>
cimg_library::CImg<T> AreaMask(cimg_library::CImg<T> origin, std::vector<Point> points)
{
	int width = origin.width(), height = origin.height();
	int x0 = width, y0 = height, x1 = 0, y1 = 0;
	cimg_library::CImg<T> mask(origin.width(), origin.height(), 1, 1, 0);
	
	//draw_line_loop
	draw_line_loop(mask, points);

	//find a point that is in the circle, we'll do this following step in order to reduce the computation
	for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++){
		x0 = x0 < (*it).x ? x0 : (*it).x;
		y0 = y0 < (*it).y ? y0 : (*it).y;
		x1 = x1 > (*it).x ? x1 : (*it).x;
		y1 = y1 > (*it).y ? y1 : (*it).y;
	}
	Point s(x0, (y1 + y0) / 2);

	int c = 0, tmp_x = points[0].x, tmp_y = points[0].y;
	while (s.x <= x1 && !(mask.atXY(s.x, s.y) == 255 && mask.atXY(s.x + 1, s.y) == 0)) {
		s.x++;
	}
	s.x++;
	
	//用栈的思路，逐渐对某点的四连通域检测并进行填充。这个方法显然8连通域会产生问题
	std::vector<Point> stack;//用了list反而变慢了一些
	stack.push_back(s);
	Point Cur = s;
	while (stack.size()){
		mask.atXY(Cur.x, Cur.y) = 255;

		if (mask.atXY(Cur.x + 1, Cur.y) == 0){
			Cur.x += 1; stack.push_back(Cur); continue;
		}
		if (mask.atXY(Cur.x, Cur.y + 1) == 0){
			Cur.y += 1; stack.push_back(Cur); continue;
		}
		if (mask.atXY(Cur.x - 1, Cur.y) == 0){
			Cur.x -= 1; stack.push_back(Cur); continue;
		}
		if (mask.atXY(Cur.x, Cur.y - 1) == 0){
			Cur.y -= 1; stack.push_back(Cur); continue;
		}
		stack.pop_back();
		if (stack.size()){
			Cur.x = stack.back().x;
			Cur.y = stack.back().y;
		}
	}

	return mask;
}//pass

//areamask with scanline
template<typename T>
cimg_library::CImg<T> _AreaMask(cimg_library::CImg<T> origin, std::vector<Point> points)
{
	cimg_library::CImg<T> mask(origin.width(), origin.height(), 1, 1, 0);
	std::vector<std::vector<float>> scanline(origin.height());//the cross points in line x;
	std::vector<Point>::iterator second = points.begin();
	std::vector<Point>::iterator first = --points.end();
	
	//generate the cross points
	for (second; second != points.end(); first = second++){ 
		if ((*second).y == (*first).y){
			scanline[(*first).y].push_back(float((*first).x));
		}
		else{
			//compute the line
			float k = float((*second).x - (*first).x) / float((*second).y - (*first).y);
			float b = float((*second).x) - k * float((*second).y);

			int del_y = (*first).y;
			while (del_y != (*second).y){
				float cross_point = float(del_y) * k + b;
				scanline[del_y].push_back(cross_point);
				if ((*second).y > (*first).y) 
					del_y++;
				else if ((*second).y < (*first).y) 
					del_y--;
			}
		}
	}
	
	//draw the mask
	for (std::vector<std::vector<float>>::size_type l = 0; l < scanline.size(); l++){
		if (scanline[l].size() > 1){
			sort(scanline[l].begin(), scanline[l].end());
			for (std::vector<std::vector<float>>::size_type it = 0; it + 1 < scanline[l].size(); it++){
				//判断线段scanline[l][it]--scanline[l][it+1]是否在多边形内部
				//由于先后两个端点在多边形上，只要计算两个点的中点是否在多边形内部即可
				if (is_inPolygon<float>((scanline[l][it] + scanline[l][it + 1]) / 2, float(l), points)) {
					for (float x = scanline[l][it]; x < scanline[l][it + 1]; x++){
						mask.atXY(int(x + 0.5), l) = 255;
					}
				}
			}
		}
	}

	return mask;
}

/**
otsu, return the computed threshold in a specific area
目前只能对灰度图做二值分割
parameter:
@origin		--original pic
@points		--ROI, 给定连续的点，按顺序连城线以确定一个封闭的环
@debug		--debug option, 0 means no information outputed
*/
int otsu(cimg_library::CImg<unsigned char> origin, std::vector<Point> points, cimg_library::CImg<unsigned char> &mask){
	int thresholdValue = 1;
	int ihist[256] = { 0 };
	int otsu_time = GetTickCount();

	int i, j, k;
	int width = origin.width(), height = origin.height();
	int n, n1, n2, gmin = 255, gmax = 0;
	double m1, m2, sum, csum, fmax, sb;
	int x0 = origin.width(), y0 = origin.height(), x1 = 0, y1 = 0;
	//缩小计算范围
	for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++){
		x0 = x0 < (*it).x ? x0 : (*it).x;
		y0 = y0 < (*it).y ? y0 : (*it).y;
		x1 = x1 > (*it).x ? x1 : (*it).x;
		y1 = y1 > (*it).y ? y1 : (*it).y;
	}

	//mask = AreaMask<unsigned char>(origin, points);
	mask = _AreaMask<unsigned char>(origin, points);

	for (i = y0; i <= y1; i++){
		for (j = x0; j < x1; j++){
			if (mask.atXY(j, i) == 255){
				ihist[origin.atXY(j, i)]++;
				if (origin.atXY(j, i) > gmax)gmax = origin.atXY(j, i);
				if (origin.atXY(j, i) < gmin)gmin = origin.atXY(j, i);
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
		if (!n1) continue;
		n2 = n - n1;
		if (n2 == 0) break;
		csum += (double)k*ihist[k];
		m1 = csum / n1;
		m2 = (sum - csum) / n2;
		sb = (double)n1*(double)n2*(m1 - m2)*(m1 - m2);
		if (sb>fmax){
			fmax = sb;
			thresholdValue = k;
		}
	}
	std::cout << "#OTSU: thresholdValue = " << thresholdValue << std::endl
		<< "gmin = " << gmin << std::endl
		<< "gmax = " << gmax << std::endl
		<< "otsu time = " << GetTickCount() - otsu_time << std::endl;

	return thresholdValue;
}//pass

/*-------------------------------------the following code find all the connected runs in a bi img----------------------------------------*/
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
			else if (origin.atXY(j - 1, i) == 255 && origin.atXY(j, i) == 0) {
				enRun.push_back(j - 1);
				masRun.push_back(enRun[NumberofRuns - 1] - stRun[NumberofRuns - 1] + 1);
			}
		}
		if (origin.atXY(origin.width() - 1, i)){
			enRun.push_back(origin.width() - 1);
			masRun.push_back(enRun[NumberofRuns - 1] - stRun[NumberofRuns - 1] + 1);
		}
	}
}

void fillRunVector(const unsigned char *origin, int width, int height, int &NumberofRuns, std::vector<int> &stRun, std::vector<int> &enRun, std::vector<int> &rowRun, std::vector<int> &masRun)
{
	for (int i = 0; i < height; i++){
		if (origin[i*width + 0] == 255){
			NumberofRuns++;
			stRun.push_back(0);
			rowRun.push_back(i);
		}
		for (int j = 1; j < width; j++){
			if (origin[i*width + (j - 1)] == 0 && origin[i*width + j] == 255){
				NumberofRuns++;
				stRun.push_back(j);
				rowRun.push_back(i);
			}
			else if (origin[i*width + (j - 1)] == 255 && origin[i*width + j] == 0) {
				enRun.push_back(j - 1);
				masRun.push_back(enRun[NumberofRuns - 1] - stRun[NumberofRuns - 1] + 1);
			}
		}
		if (origin[i*width + (width - 1)]){
			enRun.push_back(width - 1);
			masRun.push_back(enRun[NumberofRuns - 1] - stRun[NumberofRuns - 1] + 1);
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
/*------------------------------------------end-----------------------------------------------*/
//pass

//get a contour of a run
std::vector<Point> get_contour(const unsigned char *origin, int width, int height)
{
	int i = 0, j = 0;//i->y, j->x
	std::vector<Point> contour;
	int imgsize = (height + 2)*(width + 2);
	unsigned char *ori = new unsigned char[imgsize];
	memset(ori, 0, imgsize*sizeof(unsigned char));
	for (i = 0; i < height; i++){
		for (j = 0; j < width; j++){
			ori[(i + 1)*(width + 2) + (j + 1)] = origin[i*width + j];
		}
	}
	int *oriLabel = new int[imgsize];
	memset(oriLabel, 0, imgsize*sizeof(int));

	//find the start point
	for (i = 0; i < height+2; i++){
		for (j = 0; j < width + 2 && ori[i*(width + 2) + j] != 255; j++);
		if (ori[i*(width + 2) + j] == 255)
			break;
	}

	if (i == height + 2 && j == width + 2){
		std::cout << "there is no contour!" << std::endl;
		return contour;
	}

	Point start(j - 1, i - 1);
	contour.push_back(Point(j - 1, i - 1));
	oriLabel[i*(width + 2) + j] = 1;

	////设置循环计算的参数
	//Point CurP;
	//Point LastP = start;
	//int x[8] = {}
	Point CurP;
	Point LastP = start;
	int y[8] = { i, i + 1, i + 1, i + 1, i, i - 1, i - 1, i - 1 };
	int x[8] = { j + 1, j + 1, j, j - 1, j - 1, j - 1, j, j + 1 };
	int c = 7, n = 0;
	for (n = 0; n < 8; n++){
		if (ori[y[c] * (width + 2) + x[c]] == 255) {
			oriLabel[y[c] * (width + 2) + x[c]] = 1;
			break;
		}
		else
			oriLabel[y[c] * (width + 2) + x[c]] = -1;
		c = (c + 1) % 8;
	}
	if (n == 7 && ori[y[c] * (width + 2) + x[c]] == 0) {
		return contour;
	}
	else{
		CurP.y = y[c] - 1;
		CurP.x = x[c] - 1;
		contour.push_back(Point(CurP.x, CurP.y));
		//重新计算y[]，x[]，c， n
		y[0] = CurP.y + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1 + 1, y[3] = CurP.y + 1 + 1, y[4] = CurP.y + 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1 - 1, y[7] = CurP.y + 1 - 1;
		x[0] = CurP.x + 1 + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1, x[3] = CurP.x + 1 - 1, x[4] = CurP.x + 1 - 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1, x[7] = CurP.x + 1 + 1;
		n = 0, c = (c + 6) % 8;
	}
	//std::cout << start << std::endl;
	int stop = 1;
	while (stop != 0){
		for (n = 0; n < 8; n++){
			if (ori[y[c] * (width + 2) + x[c]] == 255) {
				oriLabel[y[c] * (width + 2) + x[c]] = 1;
				break;
			}
			else
				oriLabel[y[c] * (width + 2) + x[c]] = -1;
			c = (c + 1) % 8;
		}
		LastP = CurP;
		CurP.y = y[c] - 1;
		CurP.x = x[c] - 1;
		contour.push_back(Point(CurP.x, CurP.y));
		//重新计算y[]，x[]，c， n
		if (LastP.y == start.y && LastP.x == start.x){
			//for (int cc = 0; cc < 7; cc++)
			//	printf("%3d, %3d, %3d, %3d, %3d, %3d, %3d\n", oriLabel[cc * 7 + 0], oriLabel[cc * 7 + 1], oriLabel[cc * 7 + 2], oriLabel[cc * 7 + 3], oriLabel[cc * 7 + 4], oriLabel[cc * 7 + 5], oriLabel[cc * 7 + 6], oriLabel[cc * 7 + 7]);
			for (n = 0; n < 8 && (oriLabel[y[n] * (width + 2) + x[n]] != 0 || ori[y[n] * (width + 2) + x[n]] == 255); n++);
			if (n == 8 && (oriLabel[y[7] * (width + 2) + x[7]] != 0 || ori[y[7] * (width + 2) + x[7]] == 255))
				stop = 0;
			contour.pop_back();
		}
		y[0] = CurP.y + 1, y[1] = CurP.y + 1 + 1, y[2] = CurP.y + 1 + 1, y[3] = CurP.y + 1 + 1, y[4] = CurP.y + 1, y[5] = CurP.y + 1 - 1, y[6] = CurP.y + 1 - 1, y[7] = CurP.y + 1 - 1;
		x[0] = CurP.x + 1 + 1, x[1] = CurP.x + 1 + 1, x[2] = CurP.x + 1, x[3] = CurP.x + 1 - 1, x[4] = CurP.x + 1 - 1, x[5] = CurP.x + 1 - 1, x[6] = CurP.x + 1, x[7] = CurP.x + 1 + 1;
		c = (c + 6) % 8;
	}

	delete []ori;
	delete []oriLabel;
	return contour;
}//pass

//找最大的团并返回外轮廓
std::vector<Point> FindBiggestContour(cimg_library::CImg<unsigned char> origin, std::vector<Point> Points)
{
	std::vector<Point> contour;
	if (Points.size() < 3){
		std::cout << "points too less" << std::endl;
		return contour;
	}
	int width = origin.width(), height = origin.height();
	int begin = GetTickCount();
	cimg_library::CImg<unsigned char> mask;
	//compute the binary img
	int threshold = otsu(origin, Points, mask);
	unsigned char *bi_img = new unsigned char[width*height];
	for (int i = 0; i < height; i++){
		for (int j = 0; j < width; j++){
			if (origin.atXY(j, i) >= threshold && mask.atXY(j, i) == 255)
				bi_img[i*width + j] = 255;
			else
				bi_img[i*height + j] = 0;
		}
	}
	std::cout << "- generate bi pic time: " << GetTickCount() - begin << std::endl;
	begin = GetTickCount();
	//done

	//找团
	std::vector<int> stRun, enRun, rowRun, masRun;
	int NumberofRuns = 0, offset = 1, maxRun = 1;
	fillRunVector(bi_img, width, height, NumberofRuns, stRun, enRun, rowRun, masRun);
	std::vector<int> runLabels;
	std::vector<std::pair<int, int>> equivalences;
	firstPass(stRun, enRun, rowRun, NumberofRuns, runLabels, equivalences, offset);
	if (!NumberofRuns){
		std::cout << "NOTHING FOUND!!!" << std::endl
			<< "size of contour is: " << contour.size() << std::endl;
		return contour;
	}
	replaceSameLabel(runLabels, equivalences);

	int maxLabel = *max_element(runLabels.begin(), runLabels.end());
	int *MassofRuns = new int[maxLabel + 1];
	memset(MassofRuns, 0, (maxLabel + 1)*sizeof(int));

	memset(bi_img, 0, (height*width)*sizeof(unsigned char));
	for (int c = 0; c < NumberofRuns; c++){
		MassofRuns[runLabels[c]] += masRun[c];
		if (MassofRuns[runLabels[c]]>maxRun){
			maxRun = MassofRuns[runLabels[c]];
			maxLabel = runLabels[c];
		}
		int i = rowRun[c];
		for (int j = stRun[c]; j <= enRun[c]; j++){
			bi_img[i*width + j] = runLabels[c];
		}
	}
	for (int i = 0; i < height; i++)
	for (int j = 0; j < width; j++) {
		if (bi_img[i*width + j] == maxLabel)
			bi_img[i*width + j] = 255;
		else
			bi_img[i*width + j] = 0;
	}
	std::cout << "- find biggest run time: " << GetTickCount() - begin << std::endl;
	begin = GetTickCount();

	contour = get_contour(bi_img, width, height);
	std::cout << "- get contour time: " << GetTickCount() - begin << std::endl;
	begin = GetTickCount();
	
	delete []bi_img;
	delete []MassofRuns;
	return contour;
}