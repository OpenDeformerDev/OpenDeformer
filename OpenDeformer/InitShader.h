#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INITSHADER_H
#define ODER_INITSHADER_H

#include <iostream>
#include "oder.h"

namespace ODER{
	char* ReadShaderSource(const char* shaderFile){
		FILE* fp = fopen(shaderFile, "r");

		if (fp == NULL) return NULL;

		fseek(fp, 0L, SEEK_END);

		long size = ftell(fp);

		fseek(fp, 0L, SEEK_SET);
		char* buf = new char[size + 1];
		fread(buf, 1, size, fp);

		buf[size] = '\0';

		fclose(fp);

		return buf;

	}

	GLuint initShader(const char* vShaderFile, const char* fShaderFile,
		const char* gShaderFile = NULL){
		struct Shader{
			const char* filename;
			GLenum type;
			GLchar* source;
		} shaders[3] = {
				{ vShaderFile, GL_VERTEX_SHADER, NULL },
				{ fShaderFile, GL_FRAGMENT_SHADER, NULL },
				{ gShaderFile, GL_GEOMETRY_SHADER, NULL }
		};

		GLuint program = glCreateProgram();

		for (int i = 0; i < 3; i++){
			Shader& s = shaders[i];
			if (s.filename == NULL && i == 2){
				break;
			}
			s.source = ReadShaderSource(s.filename);
			if (s.source == NULL){
				std::cerr << "Fail to read " << s.filename << std::endl;
				exit(EXIT_FAILURE);
			}

			GLuint shader = glCreateShader(s.type);
			glShaderSource(shader, 1, (const GLchar**)&s.source, NULL);
			glCompileShader(shader);

			GLint compiled;
			glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
			if (!compiled){
				std::cerr << s.filename << "Fail to compile " << std::endl;
				GLint log_size;
				glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_size);
				char* log = new char[log_size];
				glGetShaderInfoLog(shader, log_size, NULL, log);
				std::cerr << log << std::endl;
				delete[] log;
				exit(EXIT_FAILURE);
			}
			delete[] s.source;

			glAttachShader(program, shader);
		}

		glLinkProgram(program);

		GLint linked;
		glGetProgramiv(program, GL_LINK_STATUS, &linked);
		if (!linked){
			std::cerr << "Fail to link" << std::endl;
			GLint log_size;
			glGetProgramiv(program, GL_INFO_LOG_LENGTH, &log_size);
			char* log = new char[log_size];
			glGetProgramInfoLog(program, log_size, NULL, log);
			std::cerr << log << std::endl;
			delete[] log;
			exit(EXIT_FAILURE);
		}

		return program;
	}
}

#endif