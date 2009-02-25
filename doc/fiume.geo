Point (1) = {0, 0, 0, 0.1};
Point (2) = {1, 1, 0, 0.1};
Point (3) = {0.2, 0.1, 0, 0.1};
Point (4) = {0.3, 0.3, 0, 0.1};
Point (5) = {0.3, 0.5, 0, 0.1};
Point (6) = {0.5, 0.6, 0, 0.1};
Point (7) = {0.7, 0.7, 0, 0.1};
Point (8) = {0.9, 0.8, 0, 0.1};
Point (9) = {1, 0.9, 0, 0.1};
Point (10) = {0.9, 1, 0, 0.1};
Point (11) = {0.8, 0.9, 0, 0.1};
Point (12) = {0.6, 0.8, 0, 0.1};
Point (13) = {0.4, 0.7, 0, 0.1};
Point (14) = {0.3, 0.6, 0, 0.1};
Point (15) = {0.2, 0.5, 0, 0.1};
Point (16) = {0.2, 0.3, 0, 0.1};
Point (17) = {0.1, 0.2, 0, 0.1};
Point (18) = {0, 0.1, 0, 0.1};
Line (1) = {1, 3};
Line (2) = {3, 4};
Line (3) = {4, 5};
Line (4) = {5, 6};
Line (5) = {6, 7};
Line (6) = {7, 8};
Line (7) = {8, 9};
Line (8) = {9, 2};
Line (9) = {2, 10};
Line (10) = {10, 11};
Line (11) = {11, 12};
Line (12) = {12, 13};
Line (13) = {13, 14};
Line (14) = {14, 15};
Line (15) = {15, 16};
Line (16) = {16, 17};
Line (17) = {17, 18};
Line (18) = {18, 1};
Line Loop(19) = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4};
Plane Surface(20) = {19};
Delete {
  Surface{20};
}
Delete {
  Line{9,10,11,12,13,14,15,16,17,1,2,3,4,5,6,7};
}
BSpline(20) = {1,3,4,5,6,7,8,9};
Delete {
  Line{20};
}
CatmullRom(20) = {1,3,4,5,6,7,8,9};
CatmullRom(21) = {2,10,11,12,13,14,15,16,17,18};
Delete {
  Line{21};
}
Delete {
  Point{10};
}
CatmullRom(21) = {2,11,12,13,14,15,16,17,18};
Translate {0,0,1} {
  Duplicata { Point{14}; }
}
Delete {
  Point{19};
}
Delete {
  Point{14};
}
Delete {
  Point{14};
}
Delete {
  Line{21};
}
Delete {
  Point{14};
}
CatmullRom(21) = {2,11,12,13,15,16,17,18};
Line Loop(22) = {20,8,21,18};
Plane Surface(23) = {22};
Delete {
  Surface{23};
}
Delete {
  Line{20};
}
Delete {
  Line{21};
}
Delete {
  Point{11,8,12,7};
}
Delete {
  Point{13,6};
}
Delete {
  Point{5,15};
}
Delete {
  Point{16,4};
}
Delete {
  Point{3,17};
}
Point(19) = {0.3,0.1,-0,0.1};
Point(20) = {0.4,0.4,-0,0.1};
Point(21) = {0.5,0.6,0,0.1};
Point(22) = {0.6,0.9,0,0.1};
Point(23) = {0.8,0.8,0,0.1};
Point(24) = {0.2,0.2,-0,0.1};
Point(25) = {0.3,0.5,0,0.1};
Point(26) = {0.4,0.7,0,0.1};
Point(27) = {0.5,1,0,0.1};
Point(28) = {0.8,0.9,0,0.1};
CatmullRom(23) = {1,19,20,21,22,23,9};
CatmullRom(24) = {18,24,25,26,27,28,2};
Line Loop(25) = {23,8,-24,18};
Plane Surface(26) = {25};
