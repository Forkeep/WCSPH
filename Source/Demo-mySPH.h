
#ifndef _DEMO_MYSPH_H_
#define _DEMO_MYSPH_H_


#ifdef PREDICTION_PCISPH
#include "PciSph.h"
#else
#ifdef IISPH
#include "IISph.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#else
#include "WcSph.h"
#endif
#endif


#ifdef PREDICTION_PCISPH
class mySPH : public PciSph {
#else
#ifdef IISPH
class mySPH : public IISph {
#else
class mySPH : public WcSph {
#endif
#endif
	glm::vec3 wheel = glm::vec3(0.0f, 7.0f, 0.0f);
	glm::vec3 wheel2 = glm::vec3(-0.0f, -7.0f, 0.0f);
	glm::vec3 wheels = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 wheel2s = glm::vec3(-0.0f, -0.0f, 0.0f);
	//��������Ҫ����������
	double rf[181][4];
	double lf[181][4];
	double rk[181][4];
	double lk[181][4];
	int flag = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;
public:

	virtual void setupScene()
	{
		if(k==1){
			
			//��ȡ�ļ����ٶ�����
			std::ifstream fin0("rf.txt", std::ios::in);
			char line[1024] = {0};
			std::string x = "";
			std::string y = "";
			std::string z = "";
			int i = 0;
			// 初始化t   数值代表从开始运行到何时开始第一次运动
			double t = 2.0f;
			//int lineflag = 0;
			while (fin0.getline(line, sizeof(line)))
			{
				cout << "����" << endl;
				//lineflag += 1;
				std::stringstream word(line);
				word >> x;
				word >> y;
				word >> z;
				/*std::cout << "x: " << x << std::endl;
			std::cout << "y: " << y << std::endl;
			std::cout << "z: " << z << std::endl;*/
				std::stringstream sa(x);
				std::stringstream sb(y);
				std::stringstream sc(z);
				double a, b, c;
				sa >> a;
				sb >> b;
				sc >> c;
				rf[i][0] = a;
				rf[i][1] = b;
				rf[i][2] = c;
				rf[i][3] = t;
				// 0.04为时间间隔，后边用到每过了这段时间，就会运动一次
				t += 0.04;
				i++;
				//std::cout << "���ǵڣ���" << lineflag << std::endl;
			}
			fin0.clear();
			fin0.close();

			//��ȡ�ļ����ٶ�����
			std::ifstream fin1("rk.txt", std::ios::in);
			/*line[1024] = { 0 };*/
			char line1[1024] = {0};
			/*std::string x = "";
		std::string y = "";
		std::string z = "";*/
			i = 0;
			t = 2.0f;
			//int lineflag = 0;
			while (fin1.getline(line1, sizeof(line1)))
			{
				cout << "����" << endl;
				//lineflag += 1;
				std::stringstream word(line1);
				word >> x;
				word >> y;
				word >> z;
				/*std::cout << "x: " << x << std::endl;
			std::cout << "y: " << y << std::endl;
			std::cout << "z: " << z << std::endl;*/
				std::stringstream sa(x);
				std::stringstream sb(y);
				std::stringstream sc(z);
				double a, b, c;
				sa >> a;
				sb >> b;
				sc >> c;
				rk[i][0] = a;
				rk[i][1] = b;
				rk[i][2] = c;
				rk[i][3] = t;
				t += 0.04;
				i++;
				//std::cout << "���ǵڣ���" << lineflag << std::endl;
			}
			fin1.clear();
			fin1.close();

			//��ȡ�ļ����ٶ�����
			std::ifstream fin2("lf.txt", std::ios::in);
			char line2[1024] = {0};
			/*std::string x = "";
		std::string y = "";
		std::string z = "";*/
			i = 0;
			t = 2.0f;
			//int lineflag = 0;
			while (fin2.getline(line2, sizeof(line2)))
			{
				cout << "����" << endl;
				//lineflag += 1;
				std::stringstream word(line2);
				word >> x;
				word >> y;
				word >> z;
				/*std::cout << "x: " << x << std::endl;
			std::cout << "y: " << y << std::endl;
			std::cout << "z: " << z << std::endl;*/
				std::stringstream sa(x);
				std::stringstream sb(y);
				std::stringstream sc(z);
				double a, b, c;
				sa >> a;
				sb >> b;
				sc >> c;
				lf[i][0] = a;
				lf[i][1] = b;
				lf[i][2] = c;
				lf[i][3] = t;
				t += 0.04;
				i++;
				//std::cout << "���ǵڣ���" << lineflag << std::endl;
			}
			fin2.clear();
			fin2.close();

			//��ȡ�ļ����ٶ�����
			std::ifstream fin3("lk.txt", std::ios::in);
			char line3[1024] = {0};
			/*std::string x = "";
		std::string y = "";
		std::string z = "";*/
			i = 0;
			t = 2.0f;
			//int lineflag = 0;
			while (fin3.getline(line3, sizeof(line3)))
			{
				cout << "����" << endl;
				//lineflag += 1;
				std::stringstream word(line3);
				word >> x;
				word >> y;
				word >> z;
				/*std::cout << "x: " << x << std::endl;
			std::cout << "y: " << y << std::endl;
			std::cout << "z: " << z << std::endl;*/
				std::stringstream sa(x);
				std::stringstream sb(y);
				std::stringstream sc(z);
				double a, b, c;
				sa >> a;
				sb >> b;
				sc >> c;
				lk[i][0] = a;
				lk[i][1] = b;
				lk[i][2] = c;
				lk[i][3] = t;
				t += 0.04;
				i++;
				//std::cout << "���ǵڣ���" << lineflag << std::endl;
			}
			fin3.clear();
			fin3.close();
		}
		


		const real_t vis = 0.05f;
		//m_TH.vis = vis;
		if (vec_t::dim == 3) {
			//����һЩϵ��//����һЩϵ��//����һЩϵ��//����һЩϵ��//����һЩϵ��//����һЩϵ��//����һЩϵ��//����һЩϵ��
			//m_TH.spacing_r = 0.05f;
			m_TH.spacing_r = 0.5f;
			//m_TH.spacing_r26 = pow(2 * m_TH.spacing_r, 6);
			m_TH.particle_volume = m_TH.spacing_r*m_TH.spacing_r*m_TH.spacing_r;
			m_TH.smoothRadius_h = 2 * m_TH.spacing_r;
			// 时间不长还和声速有关？？？？哪篇论文？
			m_TH.dt = real_t(1.0)*m_TH.spacing_r / m_TH.soundSpeed_cs * 3;
			// h是什么
			m_TH.h = 0.2f;
			m_TH.r = 1.0f;
			m_TH.bt = 1.0f;
			//m_TH.enable_vortex = 0;
			//���÷���ռ��С//���÷���ռ��С//���÷���ռ��С//���÷���ռ��С//���÷���ռ��С//���÷���ռ��С
			m_TH.spaceMin.set(-20);
			m_TH.spaceMax.set(20);
			//m_TH.spaceMax[0] = 28;
			//m_TH.spaceMax[1] = 28;

			//����2����������
			// container, solid 0
			BoundPart bp0;
			bp0.color[0] = bp0.color[1] = bp0.color[2] = 0.4f; bp0.color[3] = 1;
			bp0.position = bp0.velocity = bp0.force = vec_t::O;
			vec_t rb(3.0f);

			//vec_t container_size = vec_t(16.0f, 10.0f, 5.0f);
			//vec_t rb = container_size;
			//���ӵĴ�С�Լ����Ӻ��Ӻ��ӵĴ�С�Լ����Ӻ��Ӻ��ӵĴ�С�Լ����Ӻ��Ӻ��ӵĴ�С�Լ����Ӻ��Ӻ��ӵĴ�С�Լ����Ӻ���
			rb[0] = 13.0f;//����x
			//rb[0] = m_TH.spaceMax[0];
			//printf("%f\n",rb[0]);
			//printf("%s\n", typeid(rb).name());
			rb[1] = 4.0f;//�߶�z
			rb[2] = 8.0f;//����y
			//rb[2] = m_TH.spaceMax[2];
			//�����Ӵ�С�������棬�����߼�������
			vec_t fluidScale(1.0f);
			fluidScale[0] = rb[0];
			fluidScale[1] = rb[1] / 2;
			fluidScale[2] = rb[2];
			//����ԭ��λ��
			//addRigidCuboid(rb, bp0, 0, vis, false, glm::translate(glm::vec3(10.5f, 10.0f - m_TH.spacing_r / 2 - 2, 0)));
			//��һ����x���ڶ�����z����������y
			addRigidCuboid(rb, bp0, 0, vis, false,
				glm::translate(glm::vec3(0, -10.0f, 0)));
			//rb[0] = 0.1f;
			//rb[1] = 7.0f;
			//rb[2] = 5.0f;
			//addRigidCuboid(rb, bp0, 0, vis, false, glm::translate(glm::vec3(-1.0f, 11.0f - m_TH.spacing_r / 2 - 2, 0)));


			bp0.color[0] = bp0.color[1] = bp0.color[2] = bp0.color[3] = 0.4f;
			//bp0.color[3] = 1;
			bp0.position = bp0.velocity = bp0.force = vec_t::O;
			//������������������������������������������������������������������������������������������������������������
			rb[0] = 0.5f;//����
			rb[1] = 0.5f;//�߶�
			rb[2] = 0.5f;//����
			//��һ����x���ڶ�����z����������y
			//addRigidCuboid(rb, bp0, 0, vis, false,
			//glm::translate(glm::vec3(1.5f, 9.0f - m_TH.spacing_r / 2 - 2, 2.0f)));

			//addRigidCuboid(rb, bp0, 0, vis, false,
			//glm::translate(glm::vec3(4.7f, 9.0f - m_TH.spacing_r / 2 - 2, -2.0f)));
			//addRigidCuboid(rb, bp0, 0, vis, false,
			//glm::translate(glm::vec3(7.5f, 9.0f - m_TH.spacing_r / 2 - 2, 2.0f)));
			//��ɫ���� 1��
			addRigidCuboid(rb, bp0, 0, vis, false,
				glm::translate(glm::vec3(rf[0][0]-5.0f, rf[0][1]+10.0f, rf[0][2]-5.0f)));
			//��ɫ���� 2��
			addRigidCuboid(rb, bp0, 0, vis, false,
				glm::translate(glm::vec3(rk[0][0] - 5.0f, rk[0][1] + 10.0f, rk[0][2] - 5.0f)));
			//��ɫ���� 3��
			addRigidCuboid(rb, bp0, 0, vis, false,
				glm::translate(glm::vec3(lf[0][0] - 5.0f, lf[0][1] + 10.0f, lf[0][2] - 5.0f)));
			//��ɫ���� 4��
			addRigidCuboid(rb, bp0, 0, vis, false,
				glm::translate(glm::vec3(lk[0][0] - 5.0f, lk[0][1] + 10.0f, lk[0][2] - 5.0f)));

			// water
			FluidPart fp0;
			fp0.velocity = vec_t::O; fp0.density = 1000;
			fp0.color[1] = fp0.color[2] = 0.3f; fp0.color[0] = 0.9f; fp0.color[3] = 1;
			//���������������������������������˴�������
			vec_t lb(-0.7f), rt(0.3f); rt[1] = 2; lb[1] += 5.2f; rt[1] += 10.0f;
			//addFluidCuboid(true, 0, vec_t(-7.35f, -2.0f, -2.44f), vec_t(-3.55f, 5.5f, 2.44f), fp0, 1000, vis, 1);//��Bender 3��
			//addFluidCuboid(true, 0, vec_t(-5.2f, -2.0f, -4.44f), vec_t(-1.35f, 10.0f, 4.44f), fp0, 1000, vis, 1);
			//other set
			addFluidCuboid(true, 0, -fluidScale + vec_t(m_TH.spacing_r) - vec_t(0, fluidScale[1], 0)- vec_t(0, 10.0f, 0),
				fluidScale - vec_t(m_TH.spacing_r) - vec_t(0, 0, 0) - vec_t(0, 10.0f, 0), fp0, 1000, vis, 1);
			//addFluidCuboid(true, 0, -container_size + vec_t(m_TH.spacing_r),
				//container_size - vec_t(m_TH.spacing_r), fp0, 1000, vis, 1);
			m_TH.gravity_g[1] = real_t(-9.8);

			for (int i = 0; i < 181; i++) {
				cout << rf[i][0] << endl;
				cout << rf[i][1] << endl;
				cout << rf[i][2] << endl;
				cout << rf[i][3] << endl;
			}

		}
		else {

		}

		logInfoOfScene();
	}

protected:
	virtual void stepEvent()
	{
		float t;
		t = 1.0f;
		
			if (abs(getSystemTime()-rf[flag][3])<m_TH.dt) {
				//cout << "shijian:" << getSystemTime() << "and" << rf[flag][3] << endl;

						glm::mat4 rot;
						//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
						rot = glm::translate(glm::vec3((rf[flag+1][0]-rf[flag][0])*2, (rf[flag+1][2]-rf[flag][2])*2, (rf[flag+1][1]-rf[flag][1]))*2);
						//cout << m_TH.dt << endl;
						addSolidTransform(1, rot);
						if (flag < 179) {
								flag++;
						}
						
					}
			if (abs(getSystemTime() - rk[flag1][3]) < m_TH.dt) {
				//cout << "shijian:" << getSystemTime() << "and" << rf[flag][3] << endl;

				glm::mat4 rot;
				//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
				rot = glm::translate(glm::vec3((rk[flag1 + 1][0] - rk[flag1][0]) * 2, (rk[flag1+ 1][2] - rk[flag1][2]) * 2, (rk[flag1 + 1][1] - rk[flag1][1])) * 2);
				//cout << m_TH.dt << endl;
				addSolidTransform(2, rot);
				if (flag1 < 179) {
					flag1++;
				}

			}

			if (abs(getSystemTime() - lf[flag2][3]) < m_TH.dt) {
				//cout << "shijian:" << getSystemTime() << "and" << rf[flag][3] << endl;

				glm::mat4 rot;
				//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
				rot = glm::translate(glm::vec3((lf[flag2 + 1][0] - lf[flag2][0]) , (lf[flag2 + 1][2] - lf[flag2][2]) , (lf[flag2 + 1][1] - lf[flag2][1])) );
				//cout << m_TH.dt << endl;
				addSolidTransform(3, rot);
				if (flag2 < 179) {
					flag2++;
				}

			}

			if (abs(getSystemTime() - lf[flag3][3]) < m_TH.dt) {
				//cout << "shijian:" << getSystemTime() << "and" << rf[flag][3] << endl;

				glm::mat4 rot;
				//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
				rot = glm::translate(glm::vec3((lk[flag3 + 1][0] - lk[flag3][0]) , (lk[flag3 + 1][2] - lk[flag3][2]) , (lk[flag3 + 1][1] - lk[flag3][1])) );
				//cout << m_TH.dt << endl;
				addSolidTransform(4, rot);
				if (flag3 < 179) {
					flag3++;
				}

			}
			//cout << "a loop finish" << endl;
		
		//if (getSystemTime() < 1.0 && getSystemTime() > 0.0) {
		//	glm::mat4 rot;
		//	//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
		//	rot = glm::translate(glm::vec3(xx[1] / 10, zz[1] / 10, yy[1] / 10));
		//	//cout << m_TH.dt << endl;
		//	addSolidTransform(3, rot);
		//}
		//if (getSystemTime() < 3.0 && getSystemTime() > 1.0) {
		//	glm::mat4 rot;
		//	//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
		//	rot = glm::translate(glm::vec3(xx[5] / 10, zz[5] / 10, yy[5] / 10));
		//	//cout << m_TH.dt << endl;
		//	addSolidTransform(3, rot);
		//}

		//if (getSystemTime() < 4.0 && getSystemTime() > 3.0) {
		//	glm::mat4 rot;
		//	//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
		//	rot = glm::translate(glm::vec3(xx[10] / 10, zz[10] / 10, yy[10] / 10));
		//	//cout << m_TH.dt << endl;
		//	addSolidTransform(3, rot);
		//}
		//if (getSystemTime() < 5.0 && getSystemTime() > 4.0) {
		//	glm::mat4 rot;
		//	//rot = glm::translate(glm::vec3(xx[i], -20 * m_TH.dt, yy[i]));
		//	rot = glm::translate(glm::vec3(xx[170] / 10, zz[170] / 10, yy[170] / 10));
		//	//cout << m_TH.dt << endl;
		//	addSolidTransform(3, rot);
		//}

		
		//if (getSystemTime() < 120 && getSystemTime() > 1 * t) {
		//	glm::mat4 rot;
		//	rot = glm::translate(wheels) * glm::rotate(rot, 5 * m_TH.dt, glm::vec3(1.0f, 0.0f, 0.0f)) * glm::translate(wheel2s);
		//	addSolidTransform(3, rot);
		//}
		//if (m_TH.frameNumber % 100 == 0) {
		//	//removeFluidParts(0, 0);
		//	FluidPart fp0 = m_Fluids[0].fluidParticles[50];
		//	printf("vortex: %f, %f, %f\n", fp0.vortex[0], fp0.vortex[1], fp0.vortex[2]);
		//	printf("vortex_refinement_parameter: %f, %f, %f\n", fp0.vortex_refinement_parameter[0],
		//		fp0.vortex_refinement_parameter[1], fp0.vortex_refinement_parameter[2]);
		//	printf("vortex_velocity_refinement: %f, %f, %f\n",
		//		fp0.vortex_velocity_refinement[0], fp0.vortex_velocity_refinement[1], fp0.vortex_velocity_refinement[2]);
		//	printf("energy_loss: %f\n",
		//		fp0.energy_loss);
		//}



	}
};

#endif //#ifndef _DEMO_MYSPH_H_


