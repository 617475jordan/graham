#include <opencv2/highgui/highgui.hpp>
#include <stdio.h>  
#include <math.h>  

#include <algorithm>  
#include <vector>  
#include <opencv2/core/core.hpp>
using namespace cv;
using namespace std;

/**
凸包问题 -- Graham扫描法：
这个算法可以直接在原数据上进行运算，因此空间复杂度为O(1)。但如果将凸包的结果存储到另一数组中，
则可能在代码级别进行优化。由于在扫描凸包前要进行排序，因此时间复杂度至少为快速排序的O(nlgn)。
后面的扫描过程复杂度为O(n)，因此整个算法的复杂度为O(nlgn)。
**/

//二维点(或向量)结构体定义  
typedef vector<Point> PTARRAY;
//判断两个点(或向量)是否相等  
bool operator==(const Point &pt1, const Point &pt2) {
	return (pt1.x == pt2.x && pt1.y == pt2.y);
}
// 比较两个向量pt1和pt2分别与x轴向量(1, 0)的夹角  
bool CompareVector(const Point &pt1, const Point &pt2) {
	//求向量的模  
	float m1 = sqrt((float)(pt1.x * pt1.x + pt1.y * pt1.y));
	float m2 = sqrt((float)(pt2.x * pt2.x + pt2.y * pt2.y));
	//两个向量分别与(1, 0)求内积  
	float v1 = pt1.x / m1, v2 = pt2.x / m2;
	return (v1 > v2 || (v1 == v2 && m1 < m2));
}
//计算凸包  
void CalcConvexHull(PTARRAY &vecSrc) {
	//点集中至少应有3个点，才能构成多边形  
	if (vecSrc.size() < 3) { return; }
	//查找基点  
	Point ptBase = vecSrc.front(); //将第1个点预设为最小点  
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		//如果当前点的y值小于最小点，或y值相等，x值较小  
		if (iter->y < ptBase.y || (iter->y == ptBase.y && iter->x > ptBase.x)) {
			//将当前点作为最小点  
			ptBase = *iter;
		}
	}
	//计算出各点与基点构成的向量  
	for (PTARRAY::iterator iter = vecSrc.begin(); iter != vecSrc.end();) {
		//排除与基点相同的点，避免后面的排序计算中出现除0错误  
		if (*iter == ptBase) {
			iter = vecSrc.erase(iter);
		}
		else {
			//方向由基点到目标点  
			iter->x -= ptBase.x, iter->y -= ptBase.y;
			++iter;
		}
	}
	//按各向量与横坐标之间的夹角排序  
	sort(vecSrc.begin(), vecSrc.end(), &CompareVector);
	//删除相同的向量  
	vecSrc.erase(unique(vecSrc.begin(), vecSrc.end()), vecSrc.end());
	//计算得到首尾依次相联的向量  
	for (PTARRAY::reverse_iterator riter = vecSrc.rbegin();
		riter != vecSrc.rend() - 1; ++riter) {
		PTARRAY::reverse_iterator riNext = riter + 1;
		//向量三角形计算公式  
		riter->x -= riNext->x, riter->y -= riNext->y;
	}
	//依次删除不在凸包上的向量  
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		//回溯删除旋转方向相反的向量，使用外积判断旋转方向  
		for (PTARRAY::iterator iLast = iter - 1; iLast != vecSrc.begin();) {
			int v1 = iter->x * iLast->y, v2 = iter->y * iLast->x;
			//如果叉积小于0，则无没有逆向旋转  
			//如果叉积等于0，还需判断方向是否相逆  
			if (v1 < v2 || (v1 == v2 && iter->x * iLast->x > 0 &&
				iter->y * iLast->y > 0)) {
				break;
			}
			//删除前一个向量后，需更新当前向量，与前面的向量首尾相连  
			//向量三角形计算公式  
			iter->x += iLast->x, iter->y += iLast->y;
			iLast = (iter = vecSrc.erase(iLast)) - 1;
		}
	}
	//将所有首尾相连的向量依次累加，换算成坐标  
	vecSrc.front().x += ptBase.x, vecSrc.front().y += ptBase.y;
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		iter->x += (iter - 1)->x, iter->y += (iter - 1)->y;
	}
	//添加基点，全部的凸包计算完成  
	vecSrc.push_back(ptBase);
}

RNG rng(12345678);

int main(void)
{
	int nPtCnt = 10; //生成的随机点数  
	PTARRAY vecSrc, vecCH;
	for (int i = 0; i < nPtCnt; ++i) {
		Point ptIn = Point(rand() % 400, rand() % 400);
		vecSrc.push_back(ptIn);
		printf("%d,%d\n", ptIn.x, ptIn.y);
	}

	/// 绘出轮廓及其凸包  
	Mat drawing = Mat::zeros(500, 500, CV_8UC3);
	Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
	for (PTARRAY::iterator iter = vecSrc.begin(); iter != vecSrc.end(); ++iter)
	{
		circle(drawing, *iter, 3, color, -1, 1);
	}

	CalcConvexHull(vecSrc);
	printf("\nConvex Hull:\n");
	for (PTARRAY::iterator iter = vecSrc.begin(); iter != vecSrc.end(); ++iter) {
		printf("%d,%d\n", iter->x, iter->y);
	}

	color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
	for (PTARRAY::iterator iter = vecSrc.begin(); iter != vecSrc.end() - 1; ++iter)
	{
		line(drawing, *iter, *(iter + 1), color, 1, 8);
	}
	line(drawing, vecSrc.at(vecSrc.size() - 1), vecSrc.at(0), color, 1, 8);
	/// 把结果显示在窗体  
	namedWindow("Hull demo", CV_WINDOW_AUTOSIZE);
	imshow("Hull demo", drawing);
	waitKey(0);
	return 0;
}