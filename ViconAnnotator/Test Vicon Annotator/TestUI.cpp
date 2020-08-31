#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "GetPos.h"
#include <Windows.h>

#include "upredictor\uPredictor\ThirdParty\glew-1.9.0\include\GL\glew.h"
#include "upredictor\uPredictor\ThirdParty\glfw-3.1.1\include\GLFW\glfw3.h"
#include <gtk\gtk.h>
#include <gl\GL.h>
#include <gl\GLU.h>
using namespace std;

int main(int argc, char *argv[]) {
	string filename, F;
	GLFWwindow* window;


	GtkWidget * window2;
	gtk_init(&argc, &argv);
	window2 = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_widget_show(window2);

	GtkWidget * vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_container_add(GTK_CONTAINER(window2), vbox);

	GtkAdjustment * sliderAdjust = gtk_adjustment_new(0, 1, 64020*4.5, 1, 1, 1);
	GtkWidget * slider = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, sliderAdjust);

	gtk_box_pack_start(GTK_BOX(vbox), slider, FALSE, FALSE, 0);

	gtk_widget_show_all(window2);

	//gtk_main();
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(1024, 768, "Trajectories Animation Test", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Window Backround
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	GetPos data;
	int numFrames = 64020;
	data.parsedata2("Trajectories.csv", numFrames);

	int i = 0;
	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, 17*numFrames*3, (data.g_vertex_buffer_data - 17*numFrames*3), GL_STATIC_DRAW);

	int step = 0;
	do {

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		glPointSize(5.0f);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);
		
		glDrawArrays(GL_POINTS, gtk_adjustment_get_value(sliderAdjust), 17);
		glDisableVertexAttribArray(0);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0 && gtk_main_iteration());

	// Cleanup VBO
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	
	delete data.g_vertex_buffer_data;

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
	
	system("PAUSE");
	return 0;
}