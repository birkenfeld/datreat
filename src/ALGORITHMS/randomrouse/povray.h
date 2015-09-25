void write_povray(ap::real_2d_array phi, int N, double maxval, string name, string objn) {
	fstream file(name.c_str(), ios::out);
  file << "#declare " << objn << "=union{" << endl;
  
  int i,j;
  double size=1.0/(double)N;
  for(i=0;i<N;i++) {
	  for(j=0;j<N;j++) {
//			file << "object{Round_Box (<0,0,0>, <" << size << "," << phi(i,j)/maxval+1.0 << "," << size << ">, RBR*" << size << ", false) ";
			file << "object{box {<0,0,0>, <" << size << "," << phi(i,j)/maxval+1.0 << "," << size << ">} ";
			file << "translate<" << (double)j/(double)N-0.5 << ",-1.0," << (double)i/(double)N-0.5 << "> ";
//			file << "scale<" << size << "," << phi(i,j)/maxval << "," << size << "> ";
//			file << "translate<" << (double)j/(double)N-0.5 << ",0.0," << (double)i/(double)N-0.5 << "> ";
			//file << "texture{pigment{color rgb<" << 1.0-0.7*(double)i/(double)(N) << "," << 1.0-0.7*(double)i/(double)(N) << "," << 0.3+0.7*(double)j/(double)(N) << ">*1}}";
			file << "texture{pigment{color rgb<" << 0.0 << "," << 0.0+1.0*(double)i/(double)(N) << "," << 0.0+(double)j/(double)N << ">*1}}";
			file << "finish {ambient AMB diffuse DIFF phong PHONG reflection REFL}}";
			file << endl;
	  }
  }
  file << "}" << endl;
  file.close();
}
