#include <opencv2/highgui/highgui.hpp>
#include <stdio.h>  
#include <math.h>  

#include <algorithm>  
#include <vector>  
#include <opencv2/core/core.hpp>
using namespace cv;
using namespace std;

/**
͹������ -- Grahamɨ�跨��
����㷨����ֱ����ԭ�����Ͻ������㣬��˿ռ临�Ӷ�ΪO(1)���������͹���Ľ���洢����һ�����У�
������ڴ��뼶������Ż���������ɨ��͹��ǰҪ�����������ʱ�临�Ӷ�����Ϊ���������O(nlgn)��
�����ɨ����̸��Ӷ�ΪO(n)����������㷨�ĸ��Ӷ�ΪO(nlgn)��
**/

//��ά��(������)�ṹ�嶨��  
typedef vector<Point> PTARRAY;
//�ж�������(������)�Ƿ����  
bool operator==(const Point &pt1, const Point &pt2) {
	return (pt1.x == pt2.x && pt1.y == pt2.y);
}
// �Ƚ���������pt1��pt2�ֱ���x������(1, 0)�ļн�  
bool CompareVector(const Point &pt1, const Point &pt2) {
	//��������ģ  
	float m1 = sqrt((float)(pt1.x * pt1.x + pt1.y * pt1.y));
	float m2 = sqrt((float)(pt2.x * pt2.x + pt2.y * pt2.y));
	//���������ֱ���(1, 0)���ڻ�  
	float v1 = pt1.x / m1, v2 = pt2.x / m2;
	return (v1 > v2 || (v1 == v2 && m1 < m2));
}
//����͹��  
void CalcConvexHull(PTARRAY &vecSrc) {
	//�㼯������Ӧ��3���㣬���ܹ��ɶ����  
	if (vecSrc.size() < 3) { return; }
	//���һ���  
	Point ptBase = vecSrc.front(); //����1����Ԥ��Ϊ��С��  
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		//�����ǰ���yֵС����С�㣬��yֵ��ȣ�xֵ��С  
		if (iter->y < ptBase.y || (iter->y == ptBase.y && iter->x > ptBase.x)) {
			//����ǰ����Ϊ��С��  
			ptBase = *iter;
		}
	}
	//�������������㹹�ɵ�����  
	for (PTARRAY::iterator iter = vecSrc.begin(); iter != vecSrc.end();) {
		//�ų��������ͬ�ĵ㣬����������������г��ֳ�0����  
		if (*iter == ptBase) {
			iter = vecSrc.erase(iter);
		}
		else {
			//�����ɻ��㵽Ŀ���  
			iter->x -= ptBase.x, iter->y -= ptBase.y;
			++iter;
		}
	}
	//���������������֮��ļн�����  
	sort(vecSrc.begin(), vecSrc.end(), &CompareVector);
	//ɾ����ͬ������  
	vecSrc.erase(unique(vecSrc.begin(), vecSrc.end()), vecSrc.end());
	//����õ���β��������������  
	for (PTARRAY::reverse_iterator riter = vecSrc.rbegin();
		riter != vecSrc.rend() - 1; ++riter) {
		PTARRAY::reverse_iterator riNext = riter + 1;
		//���������μ��㹫ʽ  
		riter->x -= riNext->x, riter->y -= riNext->y;
	}
	//����ɾ������͹���ϵ�����  
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		//����ɾ����ת�����෴��������ʹ������ж���ת����  
		for (PTARRAY::iterator iLast = iter - 1; iLast != vecSrc.begin();) {
			int v1 = iter->x * iLast->y, v2 = iter->y * iLast->x;
			//������С��0������û��������ת  
			//����������0�������жϷ����Ƿ�����  
			if (v1 < v2 || (v1 == v2 && iter->x * iLast->x > 0 &&
				iter->y * iLast->y > 0)) {
				break;
			}
			//ɾ��ǰһ������������µ�ǰ��������ǰ���������β����  
			//���������μ��㹫ʽ  
			iter->x += iLast->x, iter->y += iLast->y;
			iLast = (iter = vecSrc.erase(iLast)) - 1;
		}
	}
	//��������β���������������ۼӣ����������  
	vecSrc.front().x += ptBase.x, vecSrc.front().y += ptBase.y;
	for (PTARRAY::iterator iter = vecSrc.begin() + 1; iter != vecSrc.end(); ++iter) {
		iter->x += (iter - 1)->x, iter->y += (iter - 1)->y;
	}
	//��ӻ��㣬ȫ����͹���������  
	vecSrc.push_back(ptBase);
}

RNG rng(12345678);

int main(void)
{
	int nPtCnt = 10; //���ɵ��������  
	PTARRAY vecSrc, vecCH;
	for (int i = 0; i < nPtCnt; ++i) {
		Point ptIn = Point(rand() % 400, rand() % 400);
		vecSrc.push_back(ptIn);
		printf("%d,%d\n", ptIn.x, ptIn.y);
	}

	/// �����������͹��  
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
	/// �ѽ����ʾ�ڴ���  
	namedWindow("Hull demo", CV_WINDOW_AUTOSIZE);
	imshow("Hull demo", drawing);
	waitKey(0);
	return 0;
}