// Enable printing of objects in gdb

// Print scalar 
void print(double x) {
  std::cout << x;
}

// Print vector
void print(vector<double> x) {
  std::cout << x;
}

// Print matrix
void print(matrix<double> x) {
  std::cout << x;
}

// Print array
void print(array<double> x) {
  x.print();
}

