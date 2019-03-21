/**/
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1
#include "gl_staff.h"
#include "bt_inc.h"
#include "vcg_inc.h"
#include <fstream>
#include <sstream>
#include "fUtility.h"
#include "wrap/ply/plylib.cpp"

#include "demo-mysph.h"
//#include "test.h"

mySPH sphDemo; // the global object of sph class


bool fluid_system_run=false, draw_pos=true, draw_boun=true;
bool boun_mode = true;
bool write_file = false;
bool draw_frame = true;

#define PRINT_BOOL(a) std::cout << #a": " << a << "\n"

void key_r() { fluid_system_run = !fluid_system_run; PRINT_BOOL(fluid_system_run); }
void key_b() { draw_boun=!draw_boun; PRINT_BOOL(draw_boun); }
void key_f() { draw_pos=!draw_pos; PRINT_BOOL(draw_pos); }
void key_m() { boun_mode=!boun_mode; PRINT_BOOL(boun_mode); }
void key_w() { write_file=!write_file; PRINT_BOOL(write_file); }
void key_s() { draw_frame=!draw_frame; PRINT_BOOL(draw_frame); }


void colorFunc(float reColor[4],const mySPH::FluidPart& p){
	reColor[3] = reColor[2] = 1;
	reColor[0] = reColor[1] = std::min(1.0f, p.velocity.length()/10);
}
void colorFunc_b(float reColor[4],const mySPH::BoundPart& p){
	static const float r=sphDemo.getSpacing_r(), v_max = r*r*r * 1.8f, v_min = r*r*r * 1.0f;
	reColor[3] = 1; reColor[0]=reColor[1]=reColor[2]
		= std::max(0.0f, std::min(1.0f, 1-(p.volume-v_min)/(v_max-v_min)) );
}
bool planecut(const vec_t& p){ if(p.x-p.z+0.01f>0) return true; else return false; }
void draw_unitsphere() { glCallList(1); }
bool rigid_test(int i){ return i==0; }

static const int vedio_fps = 25;
static int vedio_next_framenum = 0;

void run_and_draw()
{
	if(fluid_system_run){
		for(int i=0; i<1; ++i){
			sphDemo.runOneStep();
			if(sphDemo.getSystemTime()>=vedio_next_framenum/float(vedio_fps)) break;
		}
	}

	glMatrixMode(GL_MODELVIEW);
	if(draw_pos){
		for(size_t i=0; i<sphDemo.getNumFluids(); ++i)
			sphDemo.setFluidPartsColor(int(i), colorFunc);
		sphDemo.oglDrawFluidParts(draw_unitsphere/*, planecut*/);
	}
	
	if(draw_boun){
		float c[4] = {0,1,0.8f, 1};
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, c);
		if(boun_mode){
			//sphDemo.oglDrawSolid();
			for(size_t i=0; i<sphDemo.getNumSolids(); ++i){
				glStaff::hsl_to_rgb(float(i)/sphDemo.getNumSolids()*330+30, 1, 0.5f, c);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, c);
				sphDemo.oglDrawSolid(int(i));
			}
		}else{
			for(size_t i=0; i<sphDemo.getNumSolids(); ++i)
				sphDemo.setBoundPartsColor(int(i), colorFunc_b);
			sphDemo.oglDrawBounParts(draw_unitsphere);
		}
	}

	if(draw_frame){
		glm::vec3 s_min(sphDemo.getSpaceMin()[0], sphDemo.getSpaceMin()[1], sphDemo.getSpaceMin()[2]);
		glm::vec3 s_max(sphDemo.getSpaceMax()[0], sphDemo.getSpaceMax()[1], sphDemo.getSpaceMax()[2]);
		float c[] = {1,1,1, 1};
		wireframe(s_min, s_max, c, 1);
	}
	
	char ss[50];
	sprintf(ss,"sys time: %f",sphDemo.getSystemTime());
	glStaff::text_upperLeft(ss, 1);
	sprintf(ss,"  dt(ms): %f",1000*sphDemo.getDt());
	glStaff::text_upperLeft(ss);
	sprintf(ss,"   frame: %d",sphDemo.getFrameNumber());
	glStaff::text_upperLeft(ss);
	
#ifdef WC_TIMEADAPTIVE
	sprintf(ss,"active %%: %.2f", 100*sphDemo.wc_percentOfActive);
	glStaff::text_upperLeft(ss);
#endif

	if(write_file && fluid_system_run &&
		sphDemo.getSystemTime()>=vedio_next_framenum/float(vedio_fps) ){
		wchar_t ss[50];
		swprintf(ss, L"img/%d.png", vedio_next_framenum);
		int w, h; glStaff::get_frame_size(&w, &h);
		il_saveImgWin(ss,0,0,w,h); // png
		
		char t[50];
		for(int i=0; i<sphDemo.getNumFluids(); i++) {
			sprintf(t, "pos/%d_%d.pos", vedio_next_framenum,i);
			write_fluid_particles(t, sphDemo, i); // pos
		}
		
		glm::mat4 cube_trans;
		for(int i=0; i<sphDemo.getNumSolids(); i++) {
			sprintf(t, "mat4/%d_%d.mat4", vedio_next_framenum, i); // mat4
			sphDemo.getRigidBodyTransform(i, cube_trans);
			glStaff::save_mat_to_file(t, cube_trans);
		}

		sprintf(t, "solid-ply/%d", vedio_next_framenum);
		sphDemo.saveAllSimulatedSolidMeshToPLY(t, true, false);

		++ vedio_next_framenum;
		if(sphDemo.getSystemTime()>120.1f) exit(0);
	}else{
		if(sphDemo.getSystemTime()==0)
			vedio_next_framenum = 0;
		else
			vedio_next_framenum = int(sphDemo.getSystemTime()*vedio_fps+0.5f);
	}
}

void draw(const glm::mat4& mat_model, const glm::mat4& mat_view)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW); glLoadMatrixf(&mat_view[0][0]);
	// world
	glMultMatrixf(&mat_model[0][0]);
	//vec_t sp = (sphDemo.getSpaceMax()+sphDemo.getSpaceMin())/2;
	//glTranslatef(-sp.x, 0, -sp.z);
	run_and_draw();
}

int main(int argc, char** argv)
{
	sphDemo.setupScene();
	sphDemo.saveAllInitialSolidMeshToPLY("solid-ply/initial", true, true);
	/*
	for(size_t i=0; i<100000000L; ++i);
	for(int i=0; i<3; ++i){
	mySPH* sph = new mySPH();
	sph->runOneStep();
	sph->runOneStep();
	delete sph;
	}//*/

#ifdef PREDICTION_PCISPH
	glStaff::init_win(1000, 800, "Demo PCISPH", "C:\\Windows\\Fonts\\lucon.ttf");
#else
	glStaff::init_win(1000, 800, "Demo WCSPH", "C:\\Windows\\Fonts\\lucon.ttf");
#endif

	glStaff::init_gl();

	GLfloat vec4f[]={1, 1, 1, 1};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, vec4f); // white DIFFUSE, SPECULAR
	glLightfv(GL_LIGHT0, GL_SPECULAR, vec4f);

	glStaff::set_mat_view( glm::lookAt( glm::vec3(10,15,15), glm::vec3(0,3,0), glm::vec3(0,1,0) ) );

	glm::mat4 mat_view;
	if(glStaff::load_mat_from_file("matrix_view", mat_view))
		glStaff::set_mat_view(mat_view);
	
	glStaff::add_key_callback('R', key_r, L"run sph");
	glStaff::add_key_callback('B', key_b, L"boundary particle");
	glStaff::add_key_callback('F', key_f, L"fluid particle");
	glStaff::add_key_callback('M', key_m, L"boundary mode");
	glStaff::add_key_callback('W', key_w, L"save pos continuously every video frame");
	glStaff::add_key_callback('S', key_s, L"draw frame of spacemin and spacemax");

	glNewList(1,GL_COMPILE);
		glutSolidSphere(sphDemo.getSpacing_r()/2,10,10);
	glEndList();

	glStaff::renderLoop(draw);

	return 0;
}

#if defined(_MSC_VER) && defined(_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
class _LeakCheckStatic {
public:
	_LeakCheckStatic(){
		int tmpDbgFlag;
		tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
		tmpDbgFlag |= _CRTDBG_LEAK_CHECK_DF;
		_CrtSetDbgFlag(tmpDbgFlag);
		//_CrtSetBreakAlloc(418); // set break point at point of malloc
		char* p = new char[16]; // for test
		strcpy(p, "leak test @-@ ");
	}
} __leakObj;
#endif
