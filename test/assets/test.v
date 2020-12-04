module top (
    a , b , c , d , \result[0], \result[1] );
  input a , b , c, d;
  output \result[0], \result[1];
  wire n6, n7, n8, n9, n10;

  assign n6 = a | b;
  assign n7 = a | c;
  assign n8 = n6 | c;
  assign n9 = ~n7 & ~b;
  assign n10 = n9 | d;
  assign \result[0] = n8 & d;
  assign \result[1] = n10 & d;
endmodule
