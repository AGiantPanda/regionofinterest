#include "CImg.h"
//#include "contours.hpp"
#include "contour2.hpp"

using namespace std;
using namespace cimg_library;

#define XMAX 512
#define YMAX 512

int main(int argc, char* argv[]){

	FILE *pFile = fopen("512_512.raw", "rb");
	int imgsize = 512 * 512;
	unsigned char* pImg = new unsigned char[imgsize];
	bool readof = (imgsize == fread(pImg, sizeof(unsigned char), imgsize, pFile));
	fclose(pFile);

	CImg<unsigned char> raw(pImg, 512, 512);
	CImgDisplay disp(raw, "rawpic");
	CImgDisplay disp1(512, 512, "contour_pic", 0);
	//CImgDisplay disp2(512, 512, "original_pic", 0);
	//CImgDisplay disp3(512, 512, "scanline_pic", 0);
	//disp2.move(400, 100);
	//disp3.move(800, 100);
	//CImgDisplay bi_disp(512, 512, "bi_pic", 0);
	int mx = 0, my = 0; 
	bool redraw = false;
	bool move = false;
	vector<Point> points;
	const unsigned char white[] = { 255, 255, 255 };
	(+raw).draw_text(1, 1, "press right button to start.", white, 0, 1, 20).display(disp);

	while (!disp.is_closed() && !disp2.is_closed()){
		if (disp.mouse_x() >= 0){
			mx = disp.mouse_x(); my = disp.mouse_y();
		}
		if (disp.button() & 2){
			if (!move){
				points.clear();
				(+raw).draw_text(1, 1, "move & click, press right button to stop.\ndo not draw a circle with any cross point.", white, 0, 1, 20).display(disp);
				disp.set_key(); move = true;
			}
			else{
				(+raw).draw_text(1, 1, "press right button to start.", white, 0, 1, 20).display(disp);
				disp.set_key(); move = false; redraw = true;
			}
		}
		if ((disp.button() & 1) && move){
			points.push_back(Point(mx, my));
			if (points.size() > 1){
				CImg<unsigned char> raw_tmp = raw;
				draw_line_loop<unsigned char>(raw_tmp, points);
				(+raw_tmp).draw_text(1, 1, "move & click, press right button to stop.\ndo not draw a circle with any cross point.", white, 0, 1, 20).display(disp);
			}
		}
		if (redraw){
			int begin = GetTickCount();
			int total = GetTickCount();
			cout << "-----------begin--------------" << endl;
			vector<Point> pr = FindBiggestContour(raw, points);
			cout << "- number of points: " << points.size() << endl;
			cout << "- findbiggestcontour time: " << GetTickCount() - begin << endl;
			begin = GetTickCount();
			vector<Point>::iterator it;
			memset(pImg, 0, imgsize*sizeof(unsigned char));
			for (it = pr.begin(); it != pr.end(); it++)
				pImg[(*it).y * 512 + (*it).x] = 255;

			CImg<unsigned char> bil(pImg, 512, 512);
			cout << "- generate pic time: " << GetTickCount() - begin << endl
				<< "- totoal time: " << GetTickCount() - total << endl
				<< "------------end---------------" << endl;
			(+bil).draw_text(2, 2, "total time: %d", white, 0, 1, 13, GetTickCount() - total).display(disp1);
			//int mask_time = GetTickCount();
			//bil = AreaMask(raw, points);
			//(+bil).draw_text(2, 2, "original time: %d", white, 0, 1, 13, GetTickCount() - mask_time).display(disp2);
			//mask_time = GetTickCount();
			//bil = _AreaMask(raw, points);
			//(+bil).draw_text(2, 2, "scanline time: %d", white, 0, 1, 13, GetTickCount() - mask_time).display(disp3);
			redraw = false;
		}
		CImgDisplay::wait(disp);
	}

	return 0;
}