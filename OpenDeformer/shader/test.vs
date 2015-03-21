#version 430

in vec3 vPosition;
uniform mat4 ModalView;

void main(){
    gl_Position=ModalView*vec4(vPosition, 1.0);
}