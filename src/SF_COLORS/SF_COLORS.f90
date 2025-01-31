MODULE SF_COLORS
!SciFortran module for color definitions and manipulation
  implicit none
  !
  !RGB color type:
  type rgb_color
  !
  !A custom type for color definitions in rgb form. A large set of predefinite
  !colors are provided. Artithmetic operators are overloaded to manipulate colors
  !
     integer :: r   !Red component in the range :code:`[0:255]`
     integer :: g   !Green component in the range :code:`[0:255]`
     integer :: b   !Blue component in the range :code:`[0:255]`
  end type rgb_color
  !
  !DEFAULT RGB COLORS:
  !main
  type(rgb_color),parameter ::  black               =rgb_color(0,0,0) !:raw-html:`<span class="colored-square" style="background-color: #000000;"></span>`
  type(rgb_color),parameter ::  red                 =rgb_color(255,0,0) !:raw-html:`<span class="colored-square" style="background-color: #ff0000;"></span>`
  type(rgb_color),parameter ::  green               =rgb_color(0, 255, 0) !:raw-html:`<span class="colored-square" style="background-color: #00FF00;"></span>`
  type(rgb_color),parameter ::  orange              =rgb_color(255,193,0) !:raw-html:`<span class="colored-square" style="background-color: #FFC100;"></span>`
  type(rgb_color),parameter ::  blue                =rgb_color(0, 0, 255) !:raw-html:`<span class="colored-square" style="background-color: #0000FF;"></span>`
  type(rgb_color),parameter ::  yellow              =rgb_color(255,255,0) !:raw-html:`<span class="colored-square" style="background-color: #FFFF00;"></span>`
  type(rgb_color),parameter ::  cyan                =rgb_color(0,255,255) !:raw-html:`<span class="colored-square" style="background-color: #00FFFF;"></span>`
  type(rgb_color),parameter ::  magenta             =rgb_color(159, 0, 159) !:raw-html:`<span class="colored-square" style="background-color: #9F009F;"></span>`
  !X11/rgb.txt color codes:
  type(rgb_color),parameter :: snow                 =rgb_color(255,250,250)  !:raw-html:`<span class="colored-square" style="background-color: #FFFAFA;"></span>`
  type(rgb_color),parameter :: ghostwhite           =rgb_color(248,248,255)  !:raw-html:`<span class="colored-square" style="background-color: #F8F8FF;"></span>`
  type(rgb_color),parameter :: whitesmoke           =rgb_color(245,245,245)  !:raw-html:`<span class="colored-square" style="background-color: #F5F5F5;"></span>`
  type(rgb_color),parameter :: gainsboro            =rgb_color(220,220,220)  !:raw-html:`<span class="colored-square" style="background-color: #DCDCDC;"></span>`
  type(rgb_color),parameter :: floralwhite          =rgb_color(255,250,240)  !:raw-html:`<span class="colored-square" style="background-color: #FFFAF0;"></span>`
  type(rgb_color),parameter :: oldlace              =rgb_color(253,245,230)  !:raw-html:`<span class="colored-square" style="background-color: #FDF5E6;"></span>`
  type(rgb_color),parameter :: linen                =rgb_color(250,240,230)  !:raw-html:`<span class="colored-square" style="background-color: #FAF0E6;"></span>`
  type(rgb_color),parameter :: antiquewhite         =rgb_color(250,235,215)  !:raw-html:`<span class="colored-square" style="background-color: #FAEBD7;"></span>`
  type(rgb_color),parameter :: papayawhip           =rgb_color(255,239,213)  !:raw-html:`<span class="colored-square" style="background-color: #FFEFD5;"></span>`
  type(rgb_color),parameter :: blanchedalmond       =rgb_color(255,235,205)  !:raw-html:`<span class="colored-square" style="background-color: #FFEBCD;"></span>`
  type(rgb_color),parameter :: bisque               =rgb_color(255,228,196)  !:raw-html:`<span class="colored-square" style="background-color: #FFE4C4;"></span>`
  type(rgb_color),parameter :: peachpuff            =rgb_color(255,218,185)  !:raw-html:`<span class="colored-square" style="background-color: #FFDAB9;"></span>`
  type(rgb_color),parameter :: navajowhite          =rgb_color(255,222,173)  !:raw-html:`<span class="colored-square" style="background-color: #FFDEAD;"></span>`
  type(rgb_color),parameter :: moccasin             =rgb_color(255,228,181)  !:raw-html:`<span class="colored-square" style="background-color: #FFE4B5;"></span>`
  type(rgb_color),parameter :: cornsilk             =rgb_color(255,248,220)  !:raw-html:`<span class="colored-square" style="background-color: #FFF8DC;"></span>`
  type(rgb_color),parameter :: ivory                =rgb_color(255,255,240)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFF0;"></span>`
  type(rgb_color),parameter :: lemonchiffon         =rgb_color(255,250,205)  !:raw-html:`<span class="colored-square" style="background-color: #FFFACD;"></span>`
  type(rgb_color),parameter :: seashell             =rgb_color(255,245,238)  !:raw-html:`<span class="colored-square" style="background-color: #FFF5EE;"></span>`
  type(rgb_color),parameter :: honeydew             =rgb_color(240,255,240)  !:raw-html:`<span class="colored-square" style="background-color: #F0FFF0;"></span>`
  type(rgb_color),parameter :: mintcream            =rgb_color(245,255,250)  !:raw-html:`<span class="colored-square" style="background-color: #F5FFFA;"></span>`
  type(rgb_color),parameter :: azure                =rgb_color(240,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #F0FFFF;"></span>`
  type(rgb_color),parameter :: aliceblue            =rgb_color(240,248,255)  !:raw-html:`<span class="colored-square" style="background-color: #F0F8FF;"></span>`
  type(rgb_color),parameter :: lavender             =rgb_color(230,230,250)  !:raw-html:`<span class="colored-square" style="background-color: #E6E6FA;"></span>`
  type(rgb_color),parameter :: lavenderblush        =rgb_color(255,240,245)  !:raw-html:`<span class="colored-square" style="background-color: #FFF0F5;"></span>`
  type(rgb_color),parameter :: mistyrose            =rgb_color(255,228,225)  !:raw-html:`<span class="colored-square" style="background-color: #FFE4E1;"></span>`
  type(rgb_color),parameter :: white                =rgb_color(255,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFFF;"></span>`
  type(rgb_color),parameter :: darkslategray        =rgb_color(47,79,79)  !:raw-html:`<span class="colored-square" style="background-color: #2F4F4F;"></span>`
  type(rgb_color),parameter :: darkslategrey        =rgb_color(47,79,79)  !:raw-html:`<span class="colored-square" style="background-color: #2F4F4F;"></span>`
  type(rgb_color),parameter :: dimgray              =rgb_color(105,105,105)  !:raw-html:`<span class="colored-square" style="background-color: #696969;"></span>`
  type(rgb_color),parameter :: dimgrey              =rgb_color(105,105,105)  !:raw-html:`<span class="colored-square" style="background-color: #696969;"></span>`
  type(rgb_color),parameter :: slategray            =rgb_color(112,128,144)  !:raw-html:`<span class="colored-square" style="background-color: #708090;"></span>`
  type(rgb_color),parameter :: slategrey            =rgb_color(112,128,144)  !:raw-html:`<span class="colored-square" style="background-color: #708090;"></span>`
  type(rgb_color),parameter :: lightslategray       =rgb_color(119,136,153)  !:raw-html:`<span class="colored-square" style="background-color: #778899;"></span>`
  type(rgb_color),parameter :: lightslategrey       =rgb_color(119,136,153)  !:raw-html:`<span class="colored-square" style="background-color: #778899;"></span>`
  type(rgb_color),parameter :: gray                 =rgb_color(190,190,190)  !:raw-html:`<span class="colored-square" style="background-color: #BEBEBE;"></span>`
  type(rgb_color),parameter :: grey                 =rgb_color(190,190,190)  !:raw-html:`<span class="colored-square" style="background-color: #BEBEBE;"></span>`
  type(rgb_color),parameter :: lightgrey            =rgb_color(211,211,211)  !:raw-html:`<span class="colored-square" style="background-color: #D3D3D3;"></span>`
  type(rgb_color),parameter :: lightgray            =rgb_color(211,211,211)  !:raw-html:`<span class="colored-square" style="background-color: #D3D3D3;"></span>`
  type(rgb_color),parameter :: midnightblue         =rgb_color(25,25,112)  !:raw-html:`<span class="colored-square" style="background-color: #191970;"></span>`
  type(rgb_color),parameter :: navy                 =rgb_color(0,0,128)  !:raw-html:`<span class="colored-square" style="background-color: #000080;"></span>`
  type(rgb_color),parameter :: navyblue             =rgb_color(0,0,128)  !:raw-html:`<span class="colored-square" style="background-color: #000080;"></span>`
  type(rgb_color),parameter :: cornflowerblue       =rgb_color(100,149,237)  !:raw-html:`<span class="colored-square" style="background-color: #6495ED;"></span>`
  type(rgb_color),parameter :: darkslateblue        =rgb_color(72,61,139)  !:raw-html:`<span class="colored-square" style="background-color: #483D8B;"></span>`
  type(rgb_color),parameter :: slateblue            =rgb_color(106,90,205)  !:raw-html:`<span class="colored-square" style="background-color: #6A5ACD;"></span>`
  type(rgb_color),parameter :: mediumslateblue      =rgb_color(123,104,238)  !:raw-html:`<span class="colored-square" style="background-color: #7B68EE;"></span>`
  type(rgb_color),parameter :: lightslateblue       =rgb_color(132,112,255)  !:raw-html:`<span class="colored-square" style="background-color: #8470FF;"></span>`
  type(rgb_color),parameter :: mediumblue           =rgb_color(0,0,205)  !:raw-html:`<span class="colored-square" style="background-color: #0000CD;"></span>`
  type(rgb_color),parameter :: royalblue            =rgb_color(65,105,225)  !:raw-html:`<span class="colored-square" style="background-color: #4169E1;"></span>`
  type(rgb_color),parameter :: dodgerblue           =rgb_color(30,144,255)  !:raw-html:`<span class="colored-square" style="background-color: #1E90FF;"></span>`
  type(rgb_color),parameter :: deepskyblue          =rgb_color(0,191,255)  !:raw-html:`<span class="colored-square" style="background-color: #00BFFF;"></span>`
  type(rgb_color),parameter :: skyblue              =rgb_color(135,206,235)  !:raw-html:`<span class="colored-square" style="background-color: #87CEEB;"></span>`
  type(rgb_color),parameter :: lightskyblue         =rgb_color(135,206,250)  !:raw-html:`<span class="colored-square" style="background-color: #87CEFA;"></span>`
  type(rgb_color),parameter :: steelblue            =rgb_color(70,130,180)  !:raw-html:`<span class="colored-square" style="background-color: #4682B4;"></span>`
  type(rgb_color),parameter :: lightsteelblue       =rgb_color(176,196,222)  !:raw-html:`<span class="colored-square" style="background-color: #B0C4DE;"></span>`
  type(rgb_color),parameter :: lightblue            =rgb_color(173,216,230)  !:raw-html:`<span class="colored-square" style="background-color: #ADD8E6;"></span>`
  type(rgb_color),parameter :: powderblue           =rgb_color(176,224,230)  !:raw-html:`<span class="colored-square" style="background-color: #B0E0E6;"></span>`
  type(rgb_color),parameter :: paleturquoise        =rgb_color(175,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #AFEEEE;"></span>`
  type(rgb_color),parameter :: darkturquoise        =rgb_color(0,206,209)  !:raw-html:`<span class="colored-square" style="background-color: #00CED1;"></span>`
  type(rgb_color),parameter :: mediumturquoise      =rgb_color(72,209,204)  !:raw-html:`<span class="colored-square" style="background-color: #48D1CC;"></span>`
  type(rgb_color),parameter :: turquoise            =rgb_color(64,224,208)  !:raw-html:`<span class="colored-square" style="background-color: #40E0D0;"></span>`
  type(rgb_color),parameter :: lightcyan            =rgb_color(224,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #E0FFFF;"></span>`
  type(rgb_color),parameter :: cadetblue            =rgb_color(95,158,160)  !:raw-html:`<span class="colored-square" style="background-color: #5F9EA0;"></span>`
  type(rgb_color),parameter :: mediumaquamarine     =rgb_color(102,205,170)  !:raw-html:`<span class="colored-square" style="background-color: #66CDAA;"></span>`
  type(rgb_color),parameter :: aquamarine           =rgb_color(127,255,212)  !:raw-html:`<span class="colored-square" style="background-color: #7FFFD4;"></span>`
  type(rgb_color),parameter :: darkgreen            =rgb_color(0,100,0)  !:raw-html:`<span class="colored-square" style="background-color: #006400;"></span>`
  type(rgb_color),parameter :: darkolivegreen       =rgb_color(85,107,47)  !:raw-html:`<span class="colored-square" style="background-color: #556B2F;"></span>`
  type(rgb_color),parameter :: darkseagreen         =rgb_color(143,188,143)  !:raw-html:`<span class="colored-square" style="background-color: #8FBC8F;"></span>`
  type(rgb_color),parameter :: seagreen             =rgb_color(46,139,87)  !:raw-html:`<span class="colored-square" style="background-color: #2E8B57;"></span>`
  type(rgb_color),parameter :: mediumseagreen       =rgb_color(60,179,113)  !:raw-html:`<span class="colored-square" style="background-color: #3CB371;"></span>`
  type(rgb_color),parameter :: lightseagreen        =rgb_color(32,178,170)  !:raw-html:`<span class="colored-square" style="background-color: #20B2AA;"></span>`
  type(rgb_color),parameter :: palegreen            =rgb_color(152,251,152)  !:raw-html:`<span class="colored-square" style="background-color: #98FB98;"></span>`
  type(rgb_color),parameter :: springgreen          =rgb_color(0,255,127)  !:raw-html:`<span class="colored-square" style="background-color: #00FF7F;"></span>`
  type(rgb_color),parameter :: lawngreen            =rgb_color(124,252,0)  !:raw-html:`<span class="colored-square" style="background-color: #7CFC00;"></span>`
  type(rgb_color),parameter :: chartreuse           =rgb_color(127,255,0)  !:raw-html:`<span class="colored-square" style="background-color: #7FFF00;"></span>`
  type(rgb_color),parameter :: mediumspringgreen    =rgb_color(0,250,154)  !:raw-html:`<span class="colored-square" style="background-color: #00FA9A;"></span>`
  type(rgb_color),parameter :: greenyellow          =rgb_color(173,255,47)  !:raw-html:`<span class="colored-square" style="background-color: #ADFF2F;"></span>`
  type(rgb_color),parameter :: limegreen            =rgb_color(50,205,50)  !:raw-html:`<span class="colored-square" style="background-color: #32CD32;"></span>`
  type(rgb_color),parameter :: yellowgreen          =rgb_color(154,205,50)  !:raw-html:`<span class="colored-square" style="background-color: #9ACD32;"></span>`
  type(rgb_color),parameter :: forestgreen          =rgb_color(34,139,34)  !:raw-html:`<span class="colored-square" style="background-color: #228B22;"></span>`
  type(rgb_color),parameter :: olivedrab            =rgb_color(107,142,35)  !:raw-html:`<span class="colored-square" style="background-color: #6B8E23;"></span>`
  type(rgb_color),parameter :: darkkhaki            =rgb_color(189,183,107)  !:raw-html:`<span class="colored-square" style="background-color: #BDB76B;"></span>`
  type(rgb_color),parameter :: khaki                =rgb_color(240,230,140)  !:raw-html:`<span class="colored-square" style="background-color: #F0E68C;"></span>`
  type(rgb_color),parameter :: palegoldenrod        =rgb_color(238,232,170)  !:raw-html:`<span class="colored-square" style="background-color: #EEE8AA;"></span>`
  type(rgb_color),parameter :: lightgoldenrodyellow =rgb_color(250,250,210)  !:raw-html:`<span class="colored-square" style="background-color: #FAFAD2;"></span>`
  type(rgb_color),parameter :: lightyellow          =rgb_color(255,255,224)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFE0;"></span>`
  type(rgb_color),parameter :: gold                 =rgb_color(255,215,0)  !:raw-html:`<span class="colored-square" style="background-color: #FFD700;"></span>`
  type(rgb_color),parameter :: lightgoldenrod       =rgb_color(238,221,130)  !:raw-html:`<span class="colored-square" style="background-color: #EEDD82;"></span>`
  type(rgb_color),parameter :: goldenrod            =rgb_color(218,165,32)  !:raw-html:`<span class="colored-square" style="background-color: #DAA520;"></span>`
  type(rgb_color),parameter :: darkgoldenrod        =rgb_color(184,134,11)  !:raw-html:`<span class="colored-square" style="background-color: #B8860B;"></span>`
  type(rgb_color),parameter :: rosybrown            =rgb_color(188,143,143)  !:raw-html:`<span class="colored-square" style="background-color: #BC8F8F;"></span>`
  type(rgb_color),parameter :: indianred            =rgb_color(205,92,92)  !:raw-html:`<span class="colored-square" style="background-color: #CD5C5C;"></span>`
  type(rgb_color),parameter :: saddlebrown          =rgb_color(139,69,19)  !:raw-html:`<span class="colored-square" style="background-color: #8B4513;"></span>`
  type(rgb_color),parameter :: sienna               =rgb_color(160,82,45)  !:raw-html:`<span class="colored-square" style="background-color: #A0522D;"></span>`
  type(rgb_color),parameter :: peru                 =rgb_color(205,133,63)  !:raw-html:`<span class="colored-square" style="background-color: #CD853F;"></span>`
  type(rgb_color),parameter :: burlywood            =rgb_color(222,184,135)  !:raw-html:`<span class="colored-square" style="background-color: #DEB887;"></span>`
  type(rgb_color),parameter :: beige                =rgb_color(245,245,220)  !:raw-html:`<span class="colored-square" style="background-color: #F5F5DC;"></span>`
  type(rgb_color),parameter :: wheat                =rgb_color(245,222,179)  !:raw-html:`<span class="colored-square" style="background-color: #F5DEB3;"></span>`
  type(rgb_color),parameter :: sandybrown           =rgb_color(244,164,96)  !:raw-html:`<span class="colored-square" style="background-color: #F4A460;"></span>`
  type(rgb_color),parameter :: chocolate            =rgb_color(210,105,30)  !:raw-html:`<span class="colored-square" style="background-color: #D2691E;"></span>`
  type(rgb_color),parameter :: firebrick            =rgb_color(178,34,34)  !:raw-html:`<span class="colored-square" style="background-color: #B22222;"></span>`
  type(rgb_color),parameter :: brown                =rgb_color(165,42,42)  !:raw-html:`<span class="colored-square" style="background-color: #A52A2A;"></span>`
  type(rgb_color),parameter :: darksalmon           =rgb_color(233,150,122)  !:raw-html:`<span class="colored-square" style="background-color: #E9967A;"></span>`
  type(rgb_color),parameter :: salmon               =rgb_color(250,128,114)  !:raw-html:`<span class="colored-square" style="background-color: #FA8072;"></span>`
  type(rgb_color),parameter :: lightsalmon          =rgb_color(255,160,122)  !:raw-html:`<span class="colored-square" style="background-color: #FFA07A;"></span>`
  type(rgb_color),parameter :: darkorange           =rgb_color(255,140,0)  !:raw-html:`<span class="colored-square" style="background-color: #FF8C00;"></span>`
  type(rgb_color),parameter :: coral                =rgb_color(255,127,80)  !:raw-html:`<span class="colored-square" style="background-color: #FF7F50;"></span>`
  type(rgb_color),parameter :: lightcoral           =rgb_color(240,128,128)  !:raw-html:`<span class="colored-square" style="background-color: #F08080;"></span>`
  type(rgb_color),parameter :: tomato               =rgb_color(255,99,71)  !:raw-html:`<span class="colored-square" style="background-color: #FF6347;"></span>`
  type(rgb_color),parameter :: orangered            =rgb_color(255,69,0)  !:raw-html:`<span class="colored-square" style="background-color: #FF4500;"></span>`
  type(rgb_color),parameter :: hotpink              =rgb_color(255,105,180)  !:raw-html:`<span class="colored-square" style="background-color: #FF69B4;"></span>`
  type(rgb_color),parameter :: deeppink             =rgb_color(255,20,147)  !:raw-html:`<span class="colored-square" style="background-color: #FF1493;"></span>`
  type(rgb_color),parameter :: pink                 =rgb_color(255,192,203)  !:raw-html:`<span class="colored-square" style="background-color: #FFC0CB;"></span>`
  type(rgb_color),parameter :: lightpink            =rgb_color(255,182,193)  !:raw-html:`<span class="colored-square" style="background-color: #FFB6C1;"></span>`
  type(rgb_color),parameter :: palevioletred        =rgb_color(219,112,147)  !:raw-html:`<span class="colored-square" style="background-color: #DB7093;"></span>`
  type(rgb_color),parameter :: maroon               =rgb_color(176,48,96)  !:raw-html:`<span class="colored-square" style="background-color: #B03060;"></span>`
  type(rgb_color),parameter :: mediumvioletred      =rgb_color(199,21,133)  !:raw-html:`<span class="colored-square" style="background-color: #C71585;"></span>`
  type(rgb_color),parameter :: violetred            =rgb_color(208,32,144)  !:raw-html:`<span class="colored-square" style="background-color: #D02090;"></span>`
  type(rgb_color),parameter :: violet               =rgb_color(238,130,238)  !:raw-html:`<span class="colored-square" style="background-color: #EE82EE;"></span>`
  type(rgb_color),parameter :: plum                 =rgb_color(221,160,221)  !:raw-html:`<span class="colored-square" style="background-color: #DDA0DD;"></span>`
  type(rgb_color),parameter :: orchid               =rgb_color(218,112,214)  !:raw-html:`<span class="colored-square" style="background-color: #DA70D6;"></span>`
  type(rgb_color),parameter :: mediumorchid         =rgb_color(186,85,211)  !:raw-html:`<span class="colored-square" style="background-color: #BA55D3;"></span>`
  type(rgb_color),parameter :: darkorchid           =rgb_color(153,50,204)  !:raw-html:`<span class="colored-square" style="background-color: #9932CC;"></span>`
  type(rgb_color),parameter :: darkviolet           =rgb_color(148,0,211)  !:raw-html:`<span class="colored-square" style="background-color: #9400D3;"></span>`
  type(rgb_color),parameter :: blueviolet           =rgb_color(138,43,226)  !:raw-html:`<span class="colored-square" style="background-color: #8A2BE2;"></span>`
  type(rgb_color),parameter :: purple               =rgb_color(160,32,240)  !:raw-html:`<span class="colored-square" style="background-color: #A020F0;"></span>`
  type(rgb_color),parameter :: mediumpurple         =rgb_color(147,112,219)  !:raw-html:`<span class="colored-square" style="background-color: #9370DB;"></span>`
  type(rgb_color),parameter :: thistle              =rgb_color(216,191,216)  !:raw-html:`<span class="colored-square" style="background-color: #D8BFD8;"></span>`
  type(rgb_color),parameter :: snow1                =rgb_color(255,250,250)  !:raw-html:`<span class="colored-square" style="background-color: #FFFAFA;"></span>`
  type(rgb_color),parameter :: snow2                =rgb_color(238,233,233)  !:raw-html:`<span class="colored-square" style="background-color: #EEE9E9;"></span>`
  type(rgb_color),parameter :: snow3                =rgb_color(205,201,201)  !:raw-html:`<span class="colored-square" style="background-color: #CDC9C9;"></span>`
  type(rgb_color),parameter :: snow4                =rgb_color(139,137,137)  !:raw-html:`<span class="colored-square" style="background-color: #8B8989;"></span>`
  type(rgb_color),parameter :: seashell1            =rgb_color(255,245,238)  !:raw-html:`<span class="colored-square" style="background-color: #FFF5EE;"></span>`
  type(rgb_color),parameter :: seashell2            =rgb_color(238,229,222)  !:raw-html:`<span class="colored-square" style="background-color: #EEE5DE;"></span>`
  type(rgb_color),parameter :: seashell3            =rgb_color(205,197,191)  !:raw-html:`<span class="colored-square" style="background-color: #CDC5BF;"></span>`
  type(rgb_color),parameter :: seashell4            =rgb_color(139,134,130)  !:raw-html:`<span class="colored-square" style="background-color: #8B8682;"></span>`
  type(rgb_color),parameter :: antiquewhite1        =rgb_color(255,239,219)  !:raw-html:`<span class="colored-square" style="background-color: #FFEFDB;"></span>`
  type(rgb_color),parameter :: antiquewhite2        =rgb_color(238,223,204)  !:raw-html:`<span class="colored-square" style="background-color: #EEDFCC;"></span>`
  type(rgb_color),parameter :: antiquewhite3        =rgb_color(205,192,176)  !:raw-html:`<span class="colored-square" style="background-color: #CDC0B0;"></span>`
  type(rgb_color),parameter :: antiquewhite4        =rgb_color(139,131,120)  !:raw-html:`<span class="colored-square" style="background-color: #8B8378;"></span>`
  type(rgb_color),parameter :: bisque1              =rgb_color(255,228,196)  !:raw-html:`<span class="colored-square" style="background-color: #FFE4C4;"></span>`
  type(rgb_color),parameter :: bisque2              =rgb_color(238,213,183)  !:raw-html:`<span class="colored-square" style="background-color: #EED5B7;"></span>`
  type(rgb_color),parameter :: bisque3              =rgb_color(205,183,158)  !:raw-html:`<span class="colored-square" style="background-color: #CDB79E;"></span>`
  type(rgb_color),parameter :: bisque4              =rgb_color(139,125,107)  !:raw-html:`<span class="colored-square" style="background-color: #8B7D6B;"></span>`
  type(rgb_color),parameter :: peachpuff1           =rgb_color(255,218,185)  !:raw-html:`<span class="colored-square" style="background-color: #FFDAB9;"></span>`
  type(rgb_color),parameter :: peachpuff2           =rgb_color(238,203,173)  !:raw-html:`<span class="colored-square" style="background-color: #EECBAD;"></span>`
  type(rgb_color),parameter :: peachpuff3           =rgb_color(205,175,149)  !:raw-html:`<span class="colored-square" style="background-color: #CDAF95;"></span>`
  type(rgb_color),parameter :: peachpuff4           =rgb_color(139,119,101)  !:raw-html:`<span class="colored-square" style="background-color: #8B7765;"></span>`
  type(rgb_color),parameter :: navajowhite1         =rgb_color(255,222,173)  !:raw-html:`<span class="colored-square" style="background-color: #FFDEAD;"></span>`
  type(rgb_color),parameter :: navajowhite2         =rgb_color(238,207,161)  !:raw-html:`<span class="colored-square" style="background-color: #EECFA1;"></span>`
  type(rgb_color),parameter :: navajowhite3         =rgb_color(205,179,139)  !:raw-html:`<span class="colored-square" style="background-color: #CDB38B;"></span>`
  type(rgb_color),parameter :: navajowhite4         =rgb_color(139,121,94)  !:raw-html:`<span class="colored-square" style="background-color: #8B795E;"></span>`
  type(rgb_color),parameter :: lemonchiffon1        =rgb_color(255,250,205)  !:raw-html:`<span class="colored-square" style="background-color: #FFFACD;"></span>`
  type(rgb_color),parameter :: lemonchiffon2        =rgb_color(238,233,191)  !:raw-html:`<span class="colored-square" style="background-color: #EEE9BF;"></span>`
  type(rgb_color),parameter :: lemonchiffon3        =rgb_color(205,201,165)  !:raw-html:`<span class="colored-square" style="background-color: #CDC9A5;"></span>`
  type(rgb_color),parameter :: lemonchiffon4        =rgb_color(139,137,112)  !:raw-html:`<span class="colored-square" style="background-color: #8B8970;"></span>`
  type(rgb_color),parameter :: cornsilk1            =rgb_color(255,248,220)  !:raw-html:`<span class="colored-square" style="background-color: #FFF8DC;"></span>`
  type(rgb_color),parameter :: cornsilk2            =rgb_color(238,232,205)  !:raw-html:`<span class="colored-square" style="background-color: #EEE8CD;"></span>`
  type(rgb_color),parameter :: cornsilk3            =rgb_color(205,200,177)  !:raw-html:`<span class="colored-square" style="background-color: #CDC8B1;"></span>`
  type(rgb_color),parameter :: cornsilk4            =rgb_color(139,136,120)  !:raw-html:`<span class="colored-square" style="background-color: #8B8878;"></span>`
  type(rgb_color),parameter :: ivory1               =rgb_color(255,255,240)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFF0;"></span>`
  type(rgb_color),parameter :: ivory2               =rgb_color(238,238,224)  !:raw-html:`<span class="colored-square" style="background-color: #EEEEE0;"></span>`
  type(rgb_color),parameter :: ivory3               =rgb_color(205,205,193)  !:raw-html:`<span class="colored-square" style="background-color: #CDCDC1;"></span>`
  type(rgb_color),parameter :: ivory4               =rgb_color(139,139,131)  !:raw-html:`<span class="colored-square" style="background-color: #8B8B83;"></span>`
  type(rgb_color),parameter :: honeydew1            =rgb_color(240,255,240)  !:raw-html:`<span class="colored-square" style="background-color: #F0FFF0;"></span>`
  type(rgb_color),parameter :: honeydew2            =rgb_color(224,238,224)  !:raw-html:`<span class="colored-square" style="background-color: #E0EEE0;"></span>`
  type(rgb_color),parameter :: honeydew3            =rgb_color(193,205,193)  !:raw-html:`<span class="colored-square" style="background-color: #C1CDC1;"></span>`
  type(rgb_color),parameter :: honeydew4            =rgb_color(131,139,131)  !:raw-html:`<span class="colored-square" style="background-color: #838B83;"></span>`
  type(rgb_color),parameter :: lavenderblush1       =rgb_color(255,240,245)  !:raw-html:`<span class="colored-square" style="background-color: #FFF0F5;"></span>`
  type(rgb_color),parameter :: lavenderblush2       =rgb_color(238,224,229)  !:raw-html:`<span class="colored-square" style="background-color: #EEE0E5;"></span>`
  type(rgb_color),parameter :: lavenderblush3       =rgb_color(205,193,197)  !:raw-html:`<span class="colored-square" style="background-color: #CDC1C5;"></span>`
  type(rgb_color),parameter :: lavenderblush4       =rgb_color(139,131,134)  !:raw-html:`<span class="colored-square" style="background-color: #8B8386;"></span>`
  type(rgb_color),parameter :: mistyrose1           =rgb_color(255,228,225)  !:raw-html:`<span class="colored-square" style="background-color: #FFE4E1;"></span>`
  type(rgb_color),parameter :: mistyrose2           =rgb_color(238,213,210)  !:raw-html:`<span class="colored-square" style="background-color: #EED5D2;"></span>`
  type(rgb_color),parameter :: mistyrose3           =rgb_color(205,183,181)  !:raw-html:`<span class="colored-square" style="background-color: #CDB7B5;"></span>`
  type(rgb_color),parameter :: mistyrose4           =rgb_color(139,125,123)  !:raw-html:`<span class="colored-square" style="background-color: #8B7D7B;"></span>`
  type(rgb_color),parameter :: azure1               =rgb_color(240,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #F0FFFF;"></span>`
  type(rgb_color),parameter :: azure2               =rgb_color(224,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #E0EEEE;"></span>`
  type(rgb_color),parameter :: azure3               =rgb_color(193,205,205)  !:raw-html:`<span class="colored-square" style="background-color: #C1CDCD;"></span>`
  type(rgb_color),parameter :: azure4               =rgb_color(131,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #838B8B;"></span>`
  type(rgb_color),parameter :: slateblue1           =rgb_color(131,111,255)  !:raw-html:`<span class="colored-square" style="background-color: #836FFF;"></span>`
  type(rgb_color),parameter :: slateblue2           =rgb_color(122,103,238)  !:raw-html:`<span class="colored-square" style="background-color: #7A67EE;"></span>`
  type(rgb_color),parameter :: slateblue3           =rgb_color(105,89,205)  !:raw-html:`<span class="colored-square" style="background-color: #6959CD;"></span>`
  type(rgb_color),parameter :: slateblue4           =rgb_color(71,60,139)  !:raw-html:`<span class="colored-square" style="background-color: #473C8B;"></span>`
  type(rgb_color),parameter :: royalblue1           =rgb_color(72,118,255)  !:raw-html:`<span class="colored-square" style="background-color: #4876FF;"></span>`
  type(rgb_color),parameter :: royalblue2           =rgb_color(67,110,238)  !:raw-html:`<span class="colored-square" style="background-color: #436EEE;"></span>`
  type(rgb_color),parameter :: royalblue3           =rgb_color(58,95,205)  !:raw-html:`<span class="colored-square" style="background-color: #3A5FCD;"></span>`
  type(rgb_color),parameter :: royalblue4           =rgb_color(39,64,139)  !:raw-html:`<span class="colored-square" style="background-color: #27408B;"></span>`
  type(rgb_color),parameter :: blue1                =rgb_color(0,0,255)  !:raw-html:`<span class="colored-square" style="background-color: #0000FF;"></span>`
  type(rgb_color),parameter :: blue2                =rgb_color(0,0,238)  !:raw-html:`<span class="colored-square" style="background-color: #0000EE;"></span>`
  type(rgb_color),parameter :: blue3                =rgb_color(0,0,205)  !:raw-html:`<span class="colored-square" style="background-color: #0000CD;"></span>`
  type(rgb_color),parameter :: blue4                =rgb_color(0,0,139)  !:raw-html:`<span class="colored-square" style="background-color: #00008B;"></span>`
  type(rgb_color),parameter :: dodgerblue1          =rgb_color(30,144,255)  !:raw-html:`<span class="colored-square" style="background-color: #1E90FF;"></span>`
  type(rgb_color),parameter :: dodgerblue2          =rgb_color(28,134,238)  !:raw-html:`<span class="colored-square" style="background-color: #1C86EE;"></span>`
  type(rgb_color),parameter :: dodgerblue3          =rgb_color(24,116,205)  !:raw-html:`<span class="colored-square" style="background-color: #1874CD;"></span>`
  type(rgb_color),parameter :: dodgerblue4          =rgb_color(16,78,139)  !:raw-html:`<span class="colored-square" style="background-color: #104E8B;"></span>`
  type(rgb_color),parameter :: steelblue1           =rgb_color(99,184,255)  !:raw-html:`<span class="colored-square" style="background-color: #63B8FF;"></span>`
  type(rgb_color),parameter :: steelblue2           =rgb_color(92,172,238)  !:raw-html:`<span class="colored-square" style="background-color: #5CACEE;"></span>`
  type(rgb_color),parameter :: steelblue3           =rgb_color(79,148,205)  !:raw-html:`<span class="colored-square" style="background-color: #4F94CD;"></span>`
  type(rgb_color),parameter :: steelblue4           =rgb_color(54,100,139)  !:raw-html:`<span class="colored-square" style="background-color: #36648B;"></span>`
  type(rgb_color),parameter :: deepskyblue1         =rgb_color(0,191,255)  !:raw-html:`<span class="colored-square" style="background-color: #00BFFF;"></span>`
  type(rgb_color),parameter :: deepskyblue2         =rgb_color(0,178,238)  !:raw-html:`<span class="colored-square" style="background-color: #00B2EE;"></span>`
  type(rgb_color),parameter :: deepskyblue3         =rgb_color(0,154,205)  !:raw-html:`<span class="colored-square" style="background-color: #009ACD;"></span>`
  type(rgb_color),parameter :: deepskyblue4         =rgb_color(0,104,139)  !:raw-html:`<span class="colored-square" style="background-color: #00688B;"></span>`
  type(rgb_color),parameter :: skyblue1             =rgb_color(135,206,255)  !:raw-html:`<span class="colored-square" style="background-color: #87CEFF;"></span>`
  type(rgb_color),parameter :: skyblue2             =rgb_color(126,192,238)  !:raw-html:`<span class="colored-square" style="background-color: #7EC0EE;"></span>`
  type(rgb_color),parameter :: skyblue3             =rgb_color(108,166,205)  !:raw-html:`<span class="colored-square" style="background-color: #6CA6CD;"></span>`
  type(rgb_color),parameter :: skyblue4             =rgb_color(74,112,139)  !:raw-html:`<span class="colored-square" style="background-color: #4A708B;"></span>`
  type(rgb_color),parameter :: lightskyblue1        =rgb_color(176,226,255)  !:raw-html:`<span class="colored-square" style="background-color: #B0E2FF;"></span>`
  type(rgb_color),parameter :: lightskyblue2        =rgb_color(164,211,238)  !:raw-html:`<span class="colored-square" style="background-color: #A4D3EE;"></span>`
  type(rgb_color),parameter :: lightskyblue3        =rgb_color(141,182,205)  !:raw-html:`<span class="colored-square" style="background-color: #8DB6CD;"></span>`
  type(rgb_color),parameter :: lightskyblue4        =rgb_color(96,123,139)  !:raw-html:`<span class="colored-square" style="background-color: #607B8B;"></span>`
  type(rgb_color),parameter :: slategray1           =rgb_color(198,226,255)  !:raw-html:`<span class="colored-square" style="background-color: #C6E2FF;"></span>`
  type(rgb_color),parameter :: slategray2           =rgb_color(185,211,238)  !:raw-html:`<span class="colored-square" style="background-color: #B9D3EE;"></span>`
  type(rgb_color),parameter :: slategray3           =rgb_color(159,182,205)  !:raw-html:`<span class="colored-square" style="background-color: #9FB6CD;"></span>`
  type(rgb_color),parameter :: slategray4           =rgb_color(108,123,139)  !:raw-html:`<span class="colored-square" style="background-color: #6C7B8B;"></span>`
  type(rgb_color),parameter :: lightsteelblue1      =rgb_color(202,225,255)  !:raw-html:`<span class="colored-square" style="background-color: #CAE1FF;"></span>`
  type(rgb_color),parameter :: lightsteelblue2      =rgb_color(188,210,238)  !:raw-html:`<span class="colored-square" style="background-color: #BCD2EE;"></span>`
  type(rgb_color),parameter :: lightsteelblue3      =rgb_color(162,181,205)  !:raw-html:`<span class="colored-square" style="background-color: #A2B5CD;"></span>`
  type(rgb_color),parameter :: lightsteelblue4      =rgb_color(110,123,139)  !:raw-html:`<span class="colored-square" style="background-color: #6E7B8B;"></span>`
  type(rgb_color),parameter :: lightblue1           =rgb_color(191,239,255)  !:raw-html:`<span class="colored-square" style="background-color: #BFEFFF;"></span>`
  type(rgb_color),parameter :: lightblue2           =rgb_color(178,223,238)  !:raw-html:`<span class="colored-square" style="background-color: #B2DFEE;"></span>`
  type(rgb_color),parameter :: lightblue3           =rgb_color(154,192,205)  !:raw-html:`<span class="colored-square" style="background-color: #9AC0CD;"></span>`
  type(rgb_color),parameter :: lightblue4           =rgb_color(104,131,139)  !:raw-html:`<span class="colored-square" style="background-color: #68838B;"></span>`
  type(rgb_color),parameter :: lightcyan1           =rgb_color(224,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #E0FFFF;"></span>`
  type(rgb_color),parameter :: lightcyan2           =rgb_color(209,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #D1EEEE;"></span>`
  type(rgb_color),parameter :: lightcyan3           =rgb_color(180,205,205)  !:raw-html:`<span class="colored-square" style="background-color: #B4CDCD;"></span>`
  type(rgb_color),parameter :: lightcyan4           =rgb_color(122,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #7A8B8B;"></span>`
  type(rgb_color),parameter :: paleturquoise1       =rgb_color(187,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #BBFFFF;"></span>`
  type(rgb_color),parameter :: paleturquoise2       =rgb_color(174,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #AEEEEE;"></span>`
  type(rgb_color),parameter :: paleturquoise3       =rgb_color(150,205,205)  !:raw-html:`<span class="colored-square" style="background-color: #96CDCD;"></span>`
  type(rgb_color),parameter :: paleturquoise4       =rgb_color(102,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #668B8B;"></span>`
  type(rgb_color),parameter :: cadetblue1           =rgb_color(152,245,255)  !:raw-html:`<span class="colored-square" style="background-color: #98F5FF;"></span>`
  type(rgb_color),parameter :: cadetblue2           =rgb_color(142,229,238)  !:raw-html:`<span class="colored-square" style="background-color: #8EE5EE;"></span>`
  type(rgb_color),parameter :: cadetblue3           =rgb_color(122,197,205)  !:raw-html:`<span class="colored-square" style="background-color: #7AC5CD;"></span>`
  type(rgb_color),parameter :: cadetblue4           =rgb_color(83,134,139)  !:raw-html:`<span class="colored-square" style="background-color: #53868B;"></span>`
  type(rgb_color),parameter :: turquoise1           =rgb_color(0,245,255)  !:raw-html:`<span class="colored-square" style="background-color: #00F5FF;"></span>`
  type(rgb_color),parameter :: turquoise2           =rgb_color(0,229,238)  !:raw-html:`<span class="colored-square" style="background-color: #00E5EE;"></span>`
  type(rgb_color),parameter :: turquoise3           =rgb_color(0,197,205)  !:raw-html:`<span class="colored-square" style="background-color: #00C5CD;"></span>`
  type(rgb_color),parameter :: turquoise4           =rgb_color(0,134,139)  !:raw-html:`<span class="colored-square" style="background-color: #00868B;"></span>`
  type(rgb_color),parameter :: cyan1                =rgb_color(0,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #00FFFF;"></span>`
  type(rgb_color),parameter :: cyan2                =rgb_color(0,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #00EEEE;"></span>`
  type(rgb_color),parameter :: cyan3                =rgb_color(0,205,205)  !:raw-html:`<span class="colored-square" style="background-color: #00CDCD;"></span>`
  type(rgb_color),parameter :: cyan4                =rgb_color(0,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #008B8B;"></span>`
  type(rgb_color),parameter :: darkslategray1       =rgb_color(151,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #97FFFF;"></span>`
  type(rgb_color),parameter :: darkslategray2       =rgb_color(141,238,238)  !:raw-html:`<span class="colored-square" style="background-color: #8DEEEE;"></span>`
  type(rgb_color),parameter :: darkslategray3       =rgb_color(121,205,205)  !:raw-html:`<span class="colored-square" style="background-color: #79CDCD;"></span>`
  type(rgb_color),parameter :: darkslategray4       =rgb_color(82,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #528B8B;"></span>`
  type(rgb_color),parameter :: aquamarine1          =rgb_color(127,255,212)  !:raw-html:`<span class="colored-square" style="background-color: #7FFFD4;"></span>`
  type(rgb_color),parameter :: aquamarine2          =rgb_color(118,238,198)  !:raw-html:`<span class="colored-square" style="background-color: #76EEC6;"></span>`
  type(rgb_color),parameter :: aquamarine3          =rgb_color(102,205,170)  !:raw-html:`<span class="colored-square" style="background-color: #66CDAA;"></span>`
  type(rgb_color),parameter :: aquamarine4          =rgb_color(69,139,116)  !:raw-html:`<span class="colored-square" style="background-color: #458B74;"></span>`
  type(rgb_color),parameter :: darkseagreen1        =rgb_color(193,255,193)  !:raw-html:`<span class="colored-square" style="background-color: #C1FFC1;"></span>`
  type(rgb_color),parameter :: darkseagreen2        =rgb_color(180,238,180)  !:raw-html:`<span class="colored-square" style="background-color: #B4EEB4;"></span>`
  type(rgb_color),parameter :: darkseagreen3        =rgb_color(155,205,155)  !:raw-html:`<span class="colored-square" style="background-color: #9BCD9B;"></span>`
  type(rgb_color),parameter :: darkseagreen4        =rgb_color(105,139,105)  !:raw-html:`<span class="colored-square" style="background-color: #698B69;"></span>`
  type(rgb_color),parameter :: seagreen1            =rgb_color(84,255,159)  !:raw-html:`<span class="colored-square" style="background-color: #54FF9F;"></span>`
  type(rgb_color),parameter :: seagreen2            =rgb_color(78,238,148)  !:raw-html:`<span class="colored-square" style="background-color: #4EEE94;"></span>`
  type(rgb_color),parameter :: seagreen3            =rgb_color(67,205,128)  !:raw-html:`<span class="colored-square" style="background-color: #43CD80;"></span>`
  type(rgb_color),parameter :: seagreen4            =rgb_color(46,139,87)  !:raw-html:`<span class="colored-square" style="background-color: #2E8B57;"></span>`
  type(rgb_color),parameter :: palegreen1           =rgb_color(154,255,154)  !:raw-html:`<span class="colored-square" style="background-color: #9AFF9A;"></span>`
  type(rgb_color),parameter :: palegreen2           =rgb_color(144,238,144)  !:raw-html:`<span class="colored-square" style="background-color: #90EE90;"></span>`
  type(rgb_color),parameter :: palegreen3           =rgb_color(124,205,124)  !:raw-html:`<span class="colored-square" style="background-color: #7CCD7C;"></span>`
  type(rgb_color),parameter :: palegreen4           =rgb_color(84,139,84)  !:raw-html:`<span class="colored-square" style="background-color: #548B54;"></span>`
  type(rgb_color),parameter :: springgreen1         =rgb_color(0,255,127)  !:raw-html:`<span class="colored-square" style="background-color: #00FF7F;"></span>`
  type(rgb_color),parameter :: springgreen2         =rgb_color(0,238,118)  !:raw-html:`<span class="colored-square" style="background-color: #00EE76;"></span>`
  type(rgb_color),parameter :: springgreen3         =rgb_color(0,205,102)  !:raw-html:`<span class="colored-square" style="background-color: #00CD66;"></span>`
  type(rgb_color),parameter :: springgreen4         =rgb_color(0,139,69)  !:raw-html:`<span class="colored-square" style="background-color: #008B45;"></span>`
  type(rgb_color),parameter :: green1               =rgb_color(0,255,0)  !:raw-html:`<span class="colored-square" style="background-color: #00FF00;"></span>`
  type(rgb_color),parameter :: green2               =rgb_color(0,238,0)  !:raw-html:`<span class="colored-square" style="background-color: #00EE00;"></span>`
  type(rgb_color),parameter :: green3               =rgb_color(0,205,0)  !:raw-html:`<span class="colored-square" style="background-color: #00CD00;"></span>`
  type(rgb_color),parameter :: green4               =rgb_color(0,139,0)  !:raw-html:`<span class="colored-square" style="background-color: #008B00;"></span>`
  type(rgb_color),parameter :: chartreuse1          =rgb_color(127,255,0)  !:raw-html:`<span class="colored-square" style="background-color: #7FFF00;"></span>`
  type(rgb_color),parameter :: chartreuse2          =rgb_color(118,238,0)  !:raw-html:`<span class="colored-square" style="background-color: #76EE00;"></span>`
  type(rgb_color),parameter :: chartreuse3          =rgb_color(102,205,0)  !:raw-html:`<span class="colored-square" style="background-color: #66CD00;"></span>`
  type(rgb_color),parameter :: chartreuse4          =rgb_color(69,139,0)  !:raw-html:`<span class="colored-square" style="background-color: #458B00;"></span>`
  type(rgb_color),parameter :: olivedrab1           =rgb_color(192,255,62)  !:raw-html:`<span class="colored-square" style="background-color: #C0FF3E;"></span>`
  type(rgb_color),parameter :: olivedrab2           =rgb_color(179,238,58)  !:raw-html:`<span class="colored-square" style="background-color: #B3EE3A;"></span>`
  type(rgb_color),parameter :: olivedrab3           =rgb_color(154,205,50)  !:raw-html:`<span class="colored-square" style="background-color: #9ACD32;"></span>`
  type(rgb_color),parameter :: olivedrab4           =rgb_color(105,139,34)  !:raw-html:`<span class="colored-square" style="background-color: #698B22;"></span>`
  type(rgb_color),parameter :: darkolivegreen1      =rgb_color(202,255,112)  !:raw-html:`<span class="colored-square" style="background-color: #CAFF70;"></span>`
  type(rgb_color),parameter :: darkolivegreen2      =rgb_color(188,238,104)  !:raw-html:`<span class="colored-square" style="background-color: #BCEE68;"></span>`
  type(rgb_color),parameter :: darkolivegreen3      =rgb_color(162,205,90)  !:raw-html:`<span class="colored-square" style="background-color: #A2CD5A;"></span>`
  type(rgb_color),parameter :: darkolivegreen4      =rgb_color(110,139,61)  !:raw-html:`<span class="colored-square" style="background-color: #6E8B3D;"></span>`
  type(rgb_color),parameter :: khaki1               =rgb_color(255,246,143)  !:raw-html:`<span class="colored-square" style="background-color: #FFF68F;"></span>`
  type(rgb_color),parameter :: khaki2               =rgb_color(238,230,133)  !:raw-html:`<span class="colored-square" style="background-color: #EEE685;"></span>`
  type(rgb_color),parameter :: khaki3               =rgb_color(205,198,115)  !:raw-html:`<span class="colored-square" style="background-color: #CDC673;"></span>`
  type(rgb_color),parameter :: khaki4               =rgb_color(139,134,78)  !:raw-html:`<span class="colored-square" style="background-color: #8B864E;"></span>`
  type(rgb_color),parameter :: lightgoldenrod1      =rgb_color(255,236,139)  !:raw-html:`<span class="colored-square" style="background-color: #FFEC8B;"></span>`
  type(rgb_color),parameter :: lightgoldenrod2      =rgb_color(238,220,130)  !:raw-html:`<span class="colored-square" style="background-color: #EEDC82;"></span>`
  type(rgb_color),parameter :: lightgoldenrod3      =rgb_color(205,190,112)  !:raw-html:`<span class="colored-square" style="background-color: #CDBE70;"></span>`
  type(rgb_color),parameter :: lightgoldenrod4      =rgb_color(139,129,76)  !:raw-html:`<span class="colored-square" style="background-color: #8B814C;"></span>`
  type(rgb_color),parameter :: lightyellow1         =rgb_color(255,255,224)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFE0;"></span>`
  type(rgb_color),parameter :: lightyellow2         =rgb_color(238,238,209)  !:raw-html:`<span class="colored-square" style="background-color: #EEEED1;"></span>`
  type(rgb_color),parameter :: lightyellow3         =rgb_color(205,205,180)  !:raw-html:`<span class="colored-square" style="background-color: #CDCDB4;"></span>`
  type(rgb_color),parameter :: lightyellow4         =rgb_color(139,139,122)  !:raw-html:`<span class="colored-square" style="background-color: #8B8B7A;"></span>`
  type(rgb_color),parameter :: yellow1              =rgb_color(255,255,0)  !:raw-html:`<span class="colored-square" style="background-color: #FFFF00;"></span>`
  type(rgb_color),parameter :: yellow2              =rgb_color(238,238,0)  !:raw-html:`<span class="colored-square" style="background-color: #EEEE00;"></span>`
  type(rgb_color),parameter :: yellow3              =rgb_color(205,205,0)  !:raw-html:`<span class="colored-square" style="background-color: #CDCD00;"></span>`
  type(rgb_color),parameter :: yellow4              =rgb_color(139,139,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B8B00;"></span>`
  type(rgb_color),parameter :: gold1                =rgb_color(255,215,0)  !:raw-html:`<span class="colored-square" style="background-color: #FFD700;"></span>`
  type(rgb_color),parameter :: gold2                =rgb_color(238,201,0)  !:raw-html:`<span class="colored-square" style="background-color: #EEC900;"></span>`
  type(rgb_color),parameter :: gold3                =rgb_color(205,173,0)  !:raw-html:`<span class="colored-square" style="background-color: #CDAD00;"></span>`
  type(rgb_color),parameter :: gold4                =rgb_color(139,117,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B7500;"></span>`
  type(rgb_color),parameter :: goldenrod1           =rgb_color(255,193,37)  !:raw-html:`<span class="colored-square" style="background-color: #FFC125;"></span>`
  type(rgb_color),parameter :: goldenrod2           =rgb_color(238,180,34)  !:raw-html:`<span class="colored-square" style="background-color: #EEB422;"></span>`
  type(rgb_color),parameter :: goldenrod3           =rgb_color(205,155,29)  !:raw-html:`<span class="colored-square" style="background-color: #CD9B1D;"></span>`
  type(rgb_color),parameter :: goldenrod4           =rgb_color(139,105,20)  !:raw-html:`<span class="colored-square" style="background-color: #8B6914;"></span>`
  type(rgb_color),parameter :: darkgoldenrod1       =rgb_color(255,185,15)  !:raw-html:`<span class="colored-square" style="background-color: #FFB90F;"></span>`
  type(rgb_color),parameter :: darkgoldenrod2       =rgb_color(238,173,14)  !:raw-html:`<span class="colored-square" style="background-color: #EEAD0E;"></span>`
  type(rgb_color),parameter :: darkgoldenrod3       =rgb_color(205,149,12)  !:raw-html:`<span class="colored-square" style="background-color: #CD950C;"></span>`
  type(rgb_color),parameter :: darkgoldenrod4       =rgb_color(139,101,8)  !:raw-html:`<span class="colored-square" style="background-color: #8B6508;"></span>`
  type(rgb_color),parameter :: rosybrown1           =rgb_color(255,193,193)  !:raw-html:`<span class="colored-square" style="background-color: #FFC1C1;"></span>`
  type(rgb_color),parameter :: rosybrown2           =rgb_color(238,180,180)  !:raw-html:`<span class="colored-square" style="background-color: #EEB4B4;"></span>`
  type(rgb_color),parameter :: rosybrown3           =rgb_color(205,155,155)  !:raw-html:`<span class="colored-square" style="background-color: #CD9B9B;"></span>`
  type(rgb_color),parameter :: rosybrown4           =rgb_color(139,105,105)  !:raw-html:`<span class="colored-square" style="background-color: #8B6969;"></span>`
  type(rgb_color),parameter :: indianred1           =rgb_color(255,106,106)  !:raw-html:`<span class="colored-square" style="background-color: #FF6A6A;"></span>`
  type(rgb_color),parameter :: indianred2           =rgb_color(238,99,99)  !:raw-html:`<span class="colored-square" style="background-color: #EE6363;"></span>`
  type(rgb_color),parameter :: indianred3           =rgb_color(205,85,85)  !:raw-html:`<span class="colored-square" style="background-color: #CD5555;"></span>`
  type(rgb_color),parameter :: indianred4           =rgb_color(139,58,58)  !:raw-html:`<span class="colored-square" style="background-color: #8B3A3A;"></span>`
  type(rgb_color),parameter :: sienna1              =rgb_color(255,130,71)  !:raw-html:`<span class="colored-square" style="background-color: #FF8247;"></span>`
  type(rgb_color),parameter :: sienna2              =rgb_color(238,121,66)  !:raw-html:`<span class="colored-square" style="background-color: #EE7942;"></span>`
  type(rgb_color),parameter :: sienna3              =rgb_color(205,104,57)  !:raw-html:`<span class="colored-square" style="background-color: #CD6839;"></span>`
  type(rgb_color),parameter :: sienna4              =rgb_color(139,71,38)  !:raw-html:`<span class="colored-square" style="background-color: #8B4726;"></span>`
  type(rgb_color),parameter :: burlywood1           =rgb_color(255,211,155)  !:raw-html:`<span class="colored-square" style="background-color: #FFD39B;"></span>`
  type(rgb_color),parameter :: burlywood2           =rgb_color(238,197,145)  !:raw-html:`<span class="colored-square" style="background-color: #EEC591;"></span>`
  type(rgb_color),parameter :: burlywood3           =rgb_color(205,170,125)  !:raw-html:`<span class="colored-square" style="background-color: #CDAA7D;"></span>`
  type(rgb_color),parameter :: burlywood4           =rgb_color(139,115,85)  !:raw-html:`<span class="colored-square" style="background-color: #8B7355;"></span>`
  type(rgb_color),parameter :: wheat1               =rgb_color(255,231,186)  !:raw-html:`<span class="colored-square" style="background-color: #FFE7BA;"></span>`
  type(rgb_color),parameter :: wheat2               =rgb_color(238,216,174)  !:raw-html:`<span class="colored-square" style="background-color: #EED8AE;"></span>`
  type(rgb_color),parameter :: wheat3               =rgb_color(205,186,150)  !:raw-html:`<span class="colored-square" style="background-color: #CDBA96;"></span>`
  type(rgb_color),parameter :: wheat4               =rgb_color(139,126,102)  !:raw-html:`<span class="colored-square" style="background-color: #8B7E66;"></span>`
  type(rgb_color),parameter :: tan1                 =rgb_color(255,165,79)  !:raw-html:`<span class="colored-square" style="background-color: #FFA54F;"></span>`
  type(rgb_color),parameter :: tan2                 =rgb_color(238,154,73)  !:raw-html:`<span class="colored-square" style="background-color: #EE9A49;"></span>`
  type(rgb_color),parameter :: tan3                 =rgb_color(205,133,63)  !:raw-html:`<span class="colored-square" style="background-color: #CD853F;"></span>`
  type(rgb_color),parameter :: tan4                 =rgb_color(139,90,43)  !:raw-html:`<span class="colored-square" style="background-color: #8B5A2B;"></span>`
  type(rgb_color),parameter :: chocolate1           =rgb_color(255,127,36)  !:raw-html:`<span class="colored-square" style="background-color: #FF7F24;"></span>`
  type(rgb_color),parameter :: chocolate2           =rgb_color(238,118,33)  !:raw-html:`<span class="colored-square" style="background-color: #EE7621;"></span>`
  type(rgb_color),parameter :: chocolate3           =rgb_color(205,102,29)  !:raw-html:`<span class="colored-square" style="background-color: #CD661D;"></span>`
  type(rgb_color),parameter :: chocolate4           =rgb_color(139,69,19)  !:raw-html:`<span class="colored-square" style="background-color: #8B4513;"></span>`
  type(rgb_color),parameter :: firebrick1           =rgb_color(255,48,48)  !:raw-html:`<span class="colored-square" style="background-color: #FF3030;"></span>`
  type(rgb_color),parameter :: firebrick2           =rgb_color(238,44,44)  !:raw-html:`<span class="colored-square" style="background-color: #EE2C2C;"></span>`
  type(rgb_color),parameter :: firebrick3           =rgb_color(205,38,38)  !:raw-html:`<span class="colored-square" style="background-color: #CD2626;"></span>`
  type(rgb_color),parameter :: firebrick4           =rgb_color(139,26,26)  !:raw-html:`<span class="colored-square" style="background-color: #8B1A1A;"></span>`
  type(rgb_color),parameter :: brown1               =rgb_color(255,64,64)  !:raw-html:`<span class="colored-square" style="background-color: #FF4040;"></span>`
  type(rgb_color),parameter :: brown2               =rgb_color(238,59,59)  !:raw-html:`<span class="colored-square" style="background-color: #EE3B3B;"></span>`
  type(rgb_color),parameter :: brown3               =rgb_color(205,51,51)  !:raw-html:`<span class="colored-square" style="background-color: #CD3333;"></span>`
  type(rgb_color),parameter :: brown4               =rgb_color(139,35,35)  !:raw-html:`<span class="colored-square" style="background-color: #8B2323;"></span>`
  type(rgb_color),parameter :: salmon1              =rgb_color(255,140,105)  !:raw-html:`<span class="colored-square" style="background-color: #FF8C69;"></span>`
  type(rgb_color),parameter :: salmon2              =rgb_color(238,130,98)  !:raw-html:`<span class="colored-square" style="background-color: #EE8262;"></span>`
  type(rgb_color),parameter :: salmon3              =rgb_color(205,112,84)  !:raw-html:`<span class="colored-square" style="background-color: #CD7054;"></span>`
  type(rgb_color),parameter :: salmon4              =rgb_color(139,76,57)  !:raw-html:`<span class="colored-square" style="background-color: #8B4C39;"></span>`
  type(rgb_color),parameter :: lightsalmon1         =rgb_color(255,160,122)  !:raw-html:`<span class="colored-square" style="background-color: #FFA07A;"></span>`
  type(rgb_color),parameter :: lightsalmon2         =rgb_color(238,149,114)  !:raw-html:`<span class="colored-square" style="background-color: #EE9572;"></span>`
  type(rgb_color),parameter :: lightsalmon3         =rgb_color(205,129,98)  !:raw-html:`<span class="colored-square" style="background-color: #CD8162;"></span>`
  type(rgb_color),parameter :: lightsalmon4         =rgb_color(139,87,66)  !:raw-html:`<span class="colored-square" style="background-color: #8B5742;"></span>`
  type(rgb_color),parameter :: orange1              =rgb_color(255,165,0)  !:raw-html:`<span class="colored-square" style="background-color: #FFA500;"></span>`
  type(rgb_color),parameter :: orange2              =rgb_color(238,154,0)  !:raw-html:`<span class="colored-square" style="background-color: #EE9A00;"></span>`
  type(rgb_color),parameter :: orange3              =rgb_color(205,133,0)  !:raw-html:`<span class="colored-square" style="background-color: #CD8500;"></span>`
  type(rgb_color),parameter :: orange4              =rgb_color(139,90,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B5A00;"></span>`
  type(rgb_color),parameter :: darkorange1          =rgb_color(255,127,0)  !:raw-html:`<span class="colored-square" style="background-color: #FF7F00;"></span>`
  type(rgb_color),parameter :: darkorange2          =rgb_color(238,118,0)  !:raw-html:`<span class="colored-square" style="background-color: #EE7600;"></span>`
  type(rgb_color),parameter :: darkorange3          =rgb_color(205,102,0)  !:raw-html:`<span class="colored-square" style="background-color: #CD6600;"></span>`
  type(rgb_color),parameter :: darkorange4          =rgb_color(139,69,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B4500;"></span>`
  type(rgb_color),parameter :: coral1               =rgb_color(255,114,86)  !:raw-html:`<span class="colored-square" style="background-color: #FF7256;"></span>`
  type(rgb_color),parameter :: coral2               =rgb_color(238,106,80)  !:raw-html:`<span class="colored-square" style="background-color: #EE6A50;"></span>`
  type(rgb_color),parameter :: coral3               =rgb_color(205,91,69)  !:raw-html:`<span class="colored-square" style="background-color: #CD5B45;"></span>`
  type(rgb_color),parameter :: coral4               =rgb_color(139,62,47)  !:raw-html:`<span class="colored-square" style="background-color: #8B3E2F;"></span>`
  type(rgb_color),parameter :: tomato1              =rgb_color(255,99,71)  !:raw-html:`<span class="colored-square" style="background-color: #FF6347;"></span>`
  type(rgb_color),parameter :: tomato2              =rgb_color(238,92,66)  !:raw-html:`<span class="colored-square" style="background-color: #EE5C42;"></span>`
  type(rgb_color),parameter :: tomato3              =rgb_color(205,79,57)  !:raw-html:`<span class="colored-square" style="background-color: #CD4F39;"></span>`
  type(rgb_color),parameter :: tomato4              =rgb_color(139,54,38)  !:raw-html:`<span class="colored-square" style="background-color: #8B3626;"></span>`
  type(rgb_color),parameter :: orangered1           =rgb_color(255,69,0)  !:raw-html:`<span class="colored-square" style="background-color: #FF4500;"></span>`
  type(rgb_color),parameter :: orangered2           =rgb_color(238,64,0)  !:raw-html:`<span class="colored-square" style="background-color: #EE4000;"></span>`
  type(rgb_color),parameter :: orangered3           =rgb_color(205,55,0)  !:raw-html:`<span class="colored-square" style="background-color: #CD3700;"></span>`
  type(rgb_color),parameter :: orangered4           =rgb_color(139,37,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B2500;"></span>`
  type(rgb_color),parameter :: red1                 =rgb_color(255,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #FF0000;"></span>`
  type(rgb_color),parameter :: red2                 =rgb_color(238,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #EE0000;"></span>`
  type(rgb_color),parameter :: red3                 =rgb_color(205,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #CD0000;"></span>`
  type(rgb_color),parameter :: red4                 =rgb_color(139,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B0000;"></span>`
  type(rgb_color),parameter :: debianred            =rgb_color(215,7,81)  !:raw-html:`<span class="colored-square" style="background-color: #D70751;"></span>`
  type(rgb_color),parameter :: deeppink1            =rgb_color(255,20,147)  !:raw-html:`<span class="colored-square" style="background-color: #FF1493;"></span>`
  type(rgb_color),parameter :: deeppink2            =rgb_color(238,18,137)  !:raw-html:`<span class="colored-square" style="background-color: #EE1289;"></span>`
  type(rgb_color),parameter :: deeppink3            =rgb_color(205,16,118)  !:raw-html:`<span class="colored-square" style="background-color: #CD1076;"></span>`
  type(rgb_color),parameter :: deeppink4            =rgb_color(139,10,80)  !:raw-html:`<span class="colored-square" style="background-color: #8B0A50;"></span>`
  type(rgb_color),parameter :: hotpink1             =rgb_color(255,110,180)  !:raw-html:`<span class="colored-square" style="background-color: #FF6EB4;"></span>`
  type(rgb_color),parameter :: hotpink2             =rgb_color(238,106,167)  !:raw-html:`<span class="colored-square" style="background-color: #EE6AA7;"></span>`
  type(rgb_color),parameter :: hotpink3             =rgb_color(205,96,144)  !:raw-html:`<span class="colored-square" style="background-color: #CD6090;"></span>`
  type(rgb_color),parameter :: hotpink4             =rgb_color(139,58,98)  !:raw-html:`<span class="colored-square" style="background-color: #8B3A62;"></span>`
  type(rgb_color),parameter :: pink1                =rgb_color(255,181,197)  !:raw-html:`<span class="colored-square" style="background-color: #FFB5C5;"></span>`
  type(rgb_color),parameter :: pink2                =rgb_color(238,169,184)  !:raw-html:`<span class="colored-square" style="background-color: #EEA9B8;"></span>`
  type(rgb_color),parameter :: pink3                =rgb_color(205,145,158)  !:raw-html:`<span class="colored-square" style="background-color: #CD919E;"></span>`
  type(rgb_color),parameter :: pink4                =rgb_color(139,99,108)  !:raw-html:`<span class="colored-square" style="background-color: #8B636C;"></span>`
  type(rgb_color),parameter :: lightpink1           =rgb_color(255,174,185)  !:raw-html:`<span class="colored-square" style="background-color: #FFAEB9;"></span>`
  type(rgb_color),parameter :: lightpink2           =rgb_color(238,162,173)  !:raw-html:`<span class="colored-square" style="background-color: #EEA2AD;"></span>`
  type(rgb_color),parameter :: lightpink3           =rgb_color(205,140,149)  !:raw-html:`<span class="colored-square" style="background-color: #CD8C95;"></span>`
  type(rgb_color),parameter :: lightpink4           =rgb_color(139,95,101)  !:raw-html:`<span class="colored-square" style="background-color: #8B5F65;"></span>`
  type(rgb_color),parameter :: palevioletred1       =rgb_color(255,130,171)  !:raw-html:`<span class="colored-square" style="background-color: #FF82AB;"></span>`
  type(rgb_color),parameter :: palevioletred2       =rgb_color(238,121,159)  !:raw-html:`<span class="colored-square" style="background-color: #EE799F;"></span>`
  type(rgb_color),parameter :: palevioletred3       =rgb_color(205,104,137)  !:raw-html:`<span class="colored-square" style="background-color: #CD6889;"></span>`
  type(rgb_color),parameter :: palevioletred4       =rgb_color(139,71,93)  !:raw-html:`<span class="colored-square" style="background-color: #8B475D;"></span>`
  type(rgb_color),parameter :: maroon1              =rgb_color(255,52,179)  !:raw-html:`<span class="colored-square" style="background-color: #FF34B3;"></span>`
  type(rgb_color),parameter :: maroon2              =rgb_color(238,48,167)  !:raw-html:`<span class="colored-square" style="background-color: #EE30A7;"></span>`
  type(rgb_color),parameter :: maroon3              =rgb_color(205,41,144)  !:raw-html:`<span class="colored-square" style="background-color: #CD2990;"></span>`
  type(rgb_color),parameter :: maroon4              =rgb_color(139,28,98)  !:raw-html:`<span class="colored-square" style="background-color: #8B1C62;"></span>`
  type(rgb_color),parameter :: violetred1           =rgb_color(255,62,150)  !:raw-html:`<span class="colored-square" style="background-color: #FF3E96;"></span>`
  type(rgb_color),parameter :: violetred2           =rgb_color(238,58,140)  !:raw-html:`<span class="colored-square" style="background-color: #EE3A8C;"></span>`
  type(rgb_color),parameter :: violetred3           =rgb_color(205,50,120)  !:raw-html:`<span class="colored-square" style="background-color: #CD3278;"></span>`
  type(rgb_color),parameter :: violetred4           =rgb_color(139,34,82)  !:raw-html:`<span class="colored-square" style="background-color: #8B2252;"></span>`
  type(rgb_color),parameter :: magenta1             =rgb_color(255,0,255)  !:raw-html:`<span class="colored-square" style="background-color: #FF00FF;"></span>`
  type(rgb_color),parameter :: magenta2             =rgb_color(238,0,238)  !:raw-html:`<span class="colored-square" style="background-color: #EE00EE;"></span>`
  type(rgb_color),parameter :: magenta3             =rgb_color(205,0,205)  !:raw-html:`<span class="colored-square" style="background-color: #CD00CD;"></span>`
  type(rgb_color),parameter :: magenta4             =rgb_color(139,0,139)  !:raw-html:`<span class="colored-square" style="background-color: #8B008B;"></span>`
  type(rgb_color),parameter :: orchid1              =rgb_color(255,131,250)  !:raw-html:`<span class="colored-square" style="background-color: #FF83FA;"></span>`
  type(rgb_color),parameter :: orchid2              =rgb_color(238,122,233)  !:raw-html:`<span class="colored-square" style="background-color: #EE7AE9;"></span>`
  type(rgb_color),parameter :: orchid3              =rgb_color(205,105,201)  !:raw-html:`<span class="colored-square" style="background-color: #CD69C9;"></span>`
  type(rgb_color),parameter :: orchid4              =rgb_color(139,71,137)  !:raw-html:`<span class="colored-square" style="background-color: #8B4789;"></span>`
  type(rgb_color),parameter :: plum1                =rgb_color(255,187,255)  !:raw-html:`<span class="colored-square" style="background-color: #FFBBFF;"></span>`
  type(rgb_color),parameter :: plum2                =rgb_color(238,174,238)  !:raw-html:`<span class="colored-square" style="background-color: #EEAEEE;"></span>`
  type(rgb_color),parameter :: plum3                =rgb_color(205,150,205)  !:raw-html:`<span class="colored-square" style="background-color: #CD96CD;"></span>`
  type(rgb_color),parameter :: plum4                =rgb_color(139,102,139)  !:raw-html:`<span class="colored-square" style="background-color: #8B668B;"></span>`
  type(rgb_color),parameter :: mediumorchid1        =rgb_color(224,102,255)  !:raw-html:`<span class="colored-square" style="background-color: #E066FF;"></span>`
  type(rgb_color),parameter :: mediumorchid2        =rgb_color(209,95,238)  !:raw-html:`<span class="colored-square" style="background-color: #D15FEE;"></span>`
  type(rgb_color),parameter :: mediumorchid3        =rgb_color(180,82,205)  !:raw-html:`<span class="colored-square" style="background-color: #B452CD;"></span>`
  type(rgb_color),parameter :: mediumorchid4        =rgb_color(122,55,139)  !:raw-html:`<span class="colored-square" style="background-color: #7A378B;"></span>`
  type(rgb_color),parameter :: darkorchid1          =rgb_color(191,62,255)  !:raw-html:`<span class="colored-square" style="background-color: #BF3EFF;"></span>`
  type(rgb_color),parameter :: darkorchid2          =rgb_color(178,58,238)  !:raw-html:`<span class="colored-square" style="background-color: #B23AEE;"></span>`
  type(rgb_color),parameter :: darkorchid3          =rgb_color(154,50,205)  !:raw-html:`<span class="colored-square" style="background-color: #9A32CD;"></span>`
  type(rgb_color),parameter :: darkorchid4          =rgb_color(104,34,139)  !:raw-html:`<span class="colored-square" style="background-color: #68228B;"></span>`
  type(rgb_color),parameter :: purple1              =rgb_color(155,48,255)  !:raw-html:`<span class="colored-square" style="background-color: #9B30FF;"></span>`
  type(rgb_color),parameter :: purple2              =rgb_color(145,44,238)  !:raw-html:`<span class="colored-square" style="background-color: #912CEE;"></span>`
  type(rgb_color),parameter :: purple3              =rgb_color(125,38,205)  !:raw-html:`<span class="colored-square" style="background-color: #7D26CD;"></span>`
  type(rgb_color),parameter :: purple4              =rgb_color(85,26,139)  !:raw-html:`<span class="colored-square" style="background-color: #551A8B;"></span>`
  type(rgb_color),parameter :: mediumpurple1        =rgb_color(171,130,255)  !:raw-html:`<span class="colored-square" style="background-color: #AB82FF;"></span>`
  type(rgb_color),parameter :: mediumpurple2        =rgb_color(159,121,238)  !:raw-html:`<span class="colored-square" style="background-color: #9F79EE;"></span>`
  type(rgb_color),parameter :: mediumpurple3        =rgb_color(137,104,205)  !:raw-html:`<span class="colored-square" style="background-color: #8968CD;"></span>`
  type(rgb_color),parameter :: mediumpurple4        =rgb_color(93,71,139)  !:raw-html:`<span class="colored-square" style="background-color: #5D478B;"></span>`
  type(rgb_color),parameter :: thistle1             =rgb_color(255,225,255)  !:raw-html:`<span class="colored-square" style="background-color: #FFE1FF;"></span>`
  type(rgb_color),parameter :: thistle2             =rgb_color(238,210,238)  !:raw-html:`<span class="colored-square" style="background-color: #EED2EE;"></span>`
  type(rgb_color),parameter :: thistle3             =rgb_color(205,181,205)  !:raw-html:`<span class="colored-square" style="background-color: #CDB5CD;"></span>`
  type(rgb_color),parameter :: thistle4             =rgb_color(139,123,139)  !:raw-html:`<span class="colored-square" style="background-color: #8B7B8B;"></span>`
  type(rgb_color),parameter :: gray0                =rgb_color(0,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #000000;"></span>`
  type(rgb_color),parameter :: grey0                =rgb_color(0,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #000000;"></span>`
  type(rgb_color),parameter :: gray1                =rgb_color(3,3,3)  !:raw-html:`<span class="colored-square" style="background-color: #030303;"></span>`
  type(rgb_color),parameter :: grey1                =rgb_color(3,3,3)  !:raw-html:`<span class="colored-square" style="background-color: #030303;"></span>`
  type(rgb_color),parameter :: gray2                =rgb_color(5,5,5)  !:raw-html:`<span class="colored-square" style="background-color: #050505;"></span>`
  type(rgb_color),parameter :: grey2                =rgb_color(5,5,5)  !:raw-html:`<span class="colored-square" style="background-color: #050505;"></span>`
  type(rgb_color),parameter :: gray3                =rgb_color(8,8,8)  !:raw-html:`<span class="colored-square" style="background-color: #080808;"></span>`
  type(rgb_color),parameter :: grey3                =rgb_color(8,8,8)  !:raw-html:`<span class="colored-square" style="background-color: #080808;"></span>`
  type(rgb_color),parameter :: gray4                =rgb_color(10,10,10)  !:raw-html:`<span class="colored-square" style="background-color: #0A0A0A;"></span>`
  type(rgb_color),parameter :: grey4                =rgb_color(10,10,10)  !:raw-html:`<span class="colored-square" style="background-color: #0A0A0A;"></span>`
  type(rgb_color),parameter :: gray5                =rgb_color(13,13,13)  !:raw-html:`<span class="colored-square" style="background-color: #0D0D0D;"></span>`
  type(rgb_color),parameter :: grey5                =rgb_color(13,13,13)  !:raw-html:`<span class="colored-square" style="background-color: #0D0D0D;"></span>`
  type(rgb_color),parameter :: gray6                =rgb_color(15,15,15)  !:raw-html:`<span class="colored-square" style="background-color: #0F0F0F;"></span>`
  type(rgb_color),parameter :: grey6                =rgb_color(15,15,15)  !:raw-html:`<span class="colored-square" style="background-color: #0F0F0F;"></span>`
  type(rgb_color),parameter :: gray7                =rgb_color(18,18,18)  !:raw-html:`<span class="colored-square" style="background-color: #121212;"></span>`
  type(rgb_color),parameter :: grey7                =rgb_color(18,18,18)  !:raw-html:`<span class="colored-square" style="background-color: #121212;"></span>`
  type(rgb_color),parameter :: gray8                =rgb_color(20,20,20)  !:raw-html:`<span class="colored-square" style="background-color: #141414;"></span>`
  type(rgb_color),parameter :: grey8                =rgb_color(20,20,20)  !:raw-html:`<span class="colored-square" style="background-color: #141414;"></span>`
  type(rgb_color),parameter :: gray9                =rgb_color(23,23,23)  !:raw-html:`<span class="colored-square" style="background-color: #171717;"></span>`
  type(rgb_color),parameter :: grey9                =rgb_color(23,23,23)  !:raw-html:`<span class="colored-square" style="background-color: #171717;"></span>`
  type(rgb_color),parameter :: gray10               =rgb_color(26,26,26)  !:raw-html:`<span class="colored-square" style="background-color: #1A1A1A;"></span>`
  type(rgb_color),parameter :: grey10               =rgb_color(26,26,26)  !:raw-html:`<span class="colored-square" style="background-color: #1A1A1A;"></span>`
  type(rgb_color),parameter :: gray11               =rgb_color(28,28,28)  !:raw-html:`<span class="colored-square" style="background-color: #1C1C1C;"></span>`
  type(rgb_color),parameter :: grey11               =rgb_color(28,28,28)  !:raw-html:`<span class="colored-square" style="background-color: #1C1C1C;"></span>`
  type(rgb_color),parameter :: gray12               =rgb_color(31,31,31)  !:raw-html:`<span class="colored-square" style="background-color: #1F1F1F;"></span>`
  type(rgb_color),parameter :: grey12               =rgb_color(31,31,31)  !:raw-html:`<span class="colored-square" style="background-color: #1F1F1F;"></span>`
  type(rgb_color),parameter :: gray13               =rgb_color(33,33,33)  !:raw-html:`<span class="colored-square" style="background-color: #212121;"></span>`
  type(rgb_color),parameter :: grey13               =rgb_color(33,33,33)  !:raw-html:`<span class="colored-square" style="background-color: #212121;"></span>`
  type(rgb_color),parameter :: gray14               =rgb_color(36,36,36)  !:raw-html:`<span class="colored-square" style="background-color: #242424;"></span>`
  type(rgb_color),parameter :: grey14               =rgb_color(36,36,36)  !:raw-html:`<span class="colored-square" style="background-color: #242424;"></span>`
  type(rgb_color),parameter :: gray15               =rgb_color(38,38,38)  !:raw-html:`<span class="colored-square" style="background-color: #262626;"></span>`
  type(rgb_color),parameter :: grey15               =rgb_color(38,38,38)  !:raw-html:`<span class="colored-square" style="background-color: #262626;"></span>`
  type(rgb_color),parameter :: gray16               =rgb_color(41,41,41)  !:raw-html:`<span class="colored-square" style="background-color: #292929;"></span>`
  type(rgb_color),parameter :: grey16               =rgb_color(41,41,41)  !:raw-html:`<span class="colored-square" style="background-color: #292929;"></span>`
  type(rgb_color),parameter :: gray17               =rgb_color(43,43,43)  !:raw-html:`<span class="colored-square" style="background-color: #2B2B2B;"></span>`
  type(rgb_color),parameter :: grey17               =rgb_color(43,43,43)  !:raw-html:`<span class="colored-square" style="background-color: #2B2B2B;"></span>`
  type(rgb_color),parameter :: gray18               =rgb_color(46,46,46)  !:raw-html:`<span class="colored-square" style="background-color: #2E2E2E;"></span>`
  type(rgb_color),parameter :: grey18               =rgb_color(46,46,46)  !:raw-html:`<span class="colored-square" style="background-color: #2E2E2E;"></span>`
  type(rgb_color),parameter :: gray19               =rgb_color(48,48,48)  !:raw-html:`<span class="colored-square" style="background-color: #303030;"></span>`
  type(rgb_color),parameter :: grey19               =rgb_color(48,48,48)  !:raw-html:`<span class="colored-square" style="background-color: #303030;"></span>`
  type(rgb_color),parameter :: gray20               =rgb_color(51,51,51)  !:raw-html:`<span class="colored-square" style="background-color: #333333;"></span>`
  type(rgb_color),parameter :: grey20               =rgb_color(51,51,51)  !:raw-html:`<span class="colored-square" style="background-color: #333333;"></span>`
  type(rgb_color),parameter :: gray21               =rgb_color(54,54,54)  !:raw-html:`<span class="colored-square" style="background-color: #363636;"></span>`
  type(rgb_color),parameter :: grey21               =rgb_color(54,54,54)  !:raw-html:`<span class="colored-square" style="background-color: #363636;"></span>`
  type(rgb_color),parameter :: gray22               =rgb_color(56,56,56)  !:raw-html:`<span class="colored-square" style="background-color: #383838;"></span>`
  type(rgb_color),parameter :: grey22               =rgb_color(56,56,56)  !:raw-html:`<span class="colored-square" style="background-color: #383838;"></span>`
  type(rgb_color),parameter :: gray23               =rgb_color(59,59,59)  !:raw-html:`<span class="colored-square" style="background-color: #3B3B3B;"></span>`
  type(rgb_color),parameter :: grey23               =rgb_color(59,59,59)  !:raw-html:`<span class="colored-square" style="background-color: #3B3B3B;"></span>`
  type(rgb_color),parameter :: gray24               =rgb_color(61,61,61)  !:raw-html:`<span class="colored-square" style="background-color: #3D3D3D;"></span>`
  type(rgb_color),parameter :: grey24               =rgb_color(61,61,61)  !:raw-html:`<span class="colored-square" style="background-color: #3D3D3D;"></span>`
  type(rgb_color),parameter :: gray25               =rgb_color(64,64,64)  !:raw-html:`<span class="colored-square" style="background-color: #404040;"></span>`
  type(rgb_color),parameter :: grey25               =rgb_color(64,64,64)  !:raw-html:`<span class="colored-square" style="background-color: #404040;"></span>`
  type(rgb_color),parameter :: gray26               =rgb_color(66,66,66)  !:raw-html:`<span class="colored-square" style="background-color: #424242;"></span>`
  type(rgb_color),parameter :: grey26               =rgb_color(66,66,66)  !:raw-html:`<span class="colored-square" style="background-color: #424242;"></span>`
  type(rgb_color),parameter :: gray27               =rgb_color(69,69,69)  !:raw-html:`<span class="colored-square" style="background-color: #454545;"></span>`
  type(rgb_color),parameter :: grey27               =rgb_color(69,69,69)  !:raw-html:`<span class="colored-square" style="background-color: #454545;"></span>`
  type(rgb_color),parameter :: gray28               =rgb_color(71,71,71)  !:raw-html:`<span class="colored-square" style="background-color: #474747;"></span>`
  type(rgb_color),parameter :: grey28               =rgb_color(71,71,71)  !:raw-html:`<span class="colored-square" style="background-color: #474747;"></span>`
  type(rgb_color),parameter :: gray29               =rgb_color(74,74,74)  !:raw-html:`<span class="colored-square" style="background-color: #4A4A4A;"></span>`
  type(rgb_color),parameter :: grey29               =rgb_color(74,74,74)  !:raw-html:`<span class="colored-square" style="background-color: #4A4A4A;"></span>`
  type(rgb_color),parameter :: gray30               =rgb_color(77,77,77)  !:raw-html:`<span class="colored-square" style="background-color: #4D4D4D;"></span>`
  type(rgb_color),parameter :: grey30               =rgb_color(77,77,77)  !:raw-html:`<span class="colored-square" style="background-color: #4D4D4D;"></span>`
  type(rgb_color),parameter :: gray31               =rgb_color(79,79,79)  !:raw-html:`<span class="colored-square" style="background-color: #4F4F4F;"></span>`
  type(rgb_color),parameter :: grey31               =rgb_color(79,79,79)  !:raw-html:`<span class="colored-square" style="background-color: #4F4F4F;"></span>`
  type(rgb_color),parameter :: gray32               =rgb_color(82,82,82)  !:raw-html:`<span class="colored-square" style="background-color: #525252;"></span>`
  type(rgb_color),parameter :: grey32               =rgb_color(82,82,82)  !:raw-html:`<span class="colored-square" style="background-color: #525252;"></span>`
  type(rgb_color),parameter :: gray33               =rgb_color(84,84,84)  !:raw-html:`<span class="colored-square" style="background-color: #545454;"></span>`
  type(rgb_color),parameter :: grey33               =rgb_color(84,84,84)  !:raw-html:`<span class="colored-square" style="background-color: #545454;"></span>`
  type(rgb_color),parameter :: gray34               =rgb_color(87,87,87)  !:raw-html:`<span class="colored-square" style="background-color: #575757;"></span>`
  type(rgb_color),parameter :: grey34               =rgb_color(87,87,87)  !:raw-html:`<span class="colored-square" style="background-color: #575757;"></span>`
  type(rgb_color),parameter :: gray35               =rgb_color(89,89,89)  !:raw-html:`<span class="colored-square" style="background-color: #595959;"></span>`
  type(rgb_color),parameter :: grey35               =rgb_color(89,89,89)  !:raw-html:`<span class="colored-square" style="background-color: #595959;"></span>`
  type(rgb_color),parameter :: gray36               =rgb_color(92,92,92)  !:raw-html:`<span class="colored-square" style="background-color: #5C5C5C;"></span>`
  type(rgb_color),parameter :: grey36               =rgb_color(92,92,92)  !:raw-html:`<span class="colored-square" style="background-color: #5C5C5C;"></span>`
  type(rgb_color),parameter :: gray37               =rgb_color(94,94,94)  !:raw-html:`<span class="colored-square" style="background-color: #5E5E5E;"></span>`
  type(rgb_color),parameter :: grey37               =rgb_color(94,94,94)  !:raw-html:`<span class="colored-square" style="background-color: #5E5E5E;"></span>`
  type(rgb_color),parameter :: gray38               =rgb_color(97,97,97)  !:raw-html:`<span class="colored-square" style="background-color: #616161;"></span>`
  type(rgb_color),parameter :: grey38               =rgb_color(97,97,97)  !:raw-html:`<span class="colored-square" style="background-color: #616161;"></span>`
  type(rgb_color),parameter :: gray39               =rgb_color(99,99,99)  !:raw-html:`<span class="colored-square" style="background-color: #636363;"></span>`
  type(rgb_color),parameter :: grey39               =rgb_color(99,99,99)  !:raw-html:`<span class="colored-square" style="background-color: #636363;"></span>`
  type(rgb_color),parameter :: gray40               =rgb_color(102,102,102)  !:raw-html:`<span class="colored-square" style="background-color: #666666;"></span>`
  type(rgb_color),parameter :: grey40               =rgb_color(102,102,102)  !:raw-html:`<span class="colored-square" style="background-color: #666666;"></span>`
  type(rgb_color),parameter :: gray41               =rgb_color(105,105,105)  !:raw-html:`<span class="colored-square" style="background-color: #696969;"></span>`
  type(rgb_color),parameter :: grey41               =rgb_color(105,105,105)  !:raw-html:`<span class="colored-square" style="background-color: #696969;"></span>`
  type(rgb_color),parameter :: gray42               =rgb_color(107,107,107)  !:raw-html:`<span class="colored-square" style="background-color: #6B6B6B;"></span>`
  type(rgb_color),parameter :: grey42               =rgb_color(107,107,107)  !:raw-html:`<span class="colored-square" style="background-color: #6B6B6B;"></span>`
  type(rgb_color),parameter :: gray43               =rgb_color(110,110,110)  !:raw-html:`<span class="colored-square" style="background-color: #6E6E6E;"></span>`
  type(rgb_color),parameter :: grey43               =rgb_color(110,110,110)  !:raw-html:`<span class="colored-square" style="background-color: #6E6E6E;"></span>`
  type(rgb_color),parameter :: gray44               =rgb_color(112,112,112)  !:raw-html:`<span class="colored-square" style="background-color: #707070;"></span>`
  type(rgb_color),parameter :: grey44               =rgb_color(112,112,112)  !:raw-html:`<span class="colored-square" style="background-color: #707070;"></span>`
  type(rgb_color),parameter :: gray45               =rgb_color(115,115,115)  !:raw-html:`<span class="colored-square" style="background-color: #737373;"></span>`
  type(rgb_color),parameter :: grey45               =rgb_color(115,115,115)  !:raw-html:`<span class="colored-square" style="background-color: #737373;"></span>`
  type(rgb_color),parameter :: gray46               =rgb_color(117,117,117)  !:raw-html:`<span class="colored-square" style="background-color: #757575;"></span>`
  type(rgb_color),parameter :: grey46               =rgb_color(117,117,117)  !:raw-html:`<span class="colored-square" style="background-color: #757575;"></span>`
  type(rgb_color),parameter :: gray47               =rgb_color(120,120,120)  !:raw-html:`<span class="colored-square" style="background-color: #787878;"></span>`
  type(rgb_color),parameter :: grey47               =rgb_color(120,120,120)  !:raw-html:`<span class="colored-square" style="background-color: #787878;"></span>`
  type(rgb_color),parameter :: gray48               =rgb_color(122,122,122)  !:raw-html:`<span class="colored-square" style="background-color: #7A7A7A;"></span>`
  type(rgb_color),parameter :: grey48               =rgb_color(122,122,122)  !:raw-html:`<span class="colored-square" style="background-color: #7A7A7A;"></span>`
  type(rgb_color),parameter :: gray49               =rgb_color(125,125,125)  !:raw-html:`<span class="colored-square" style="background-color: #7D7D7D;"></span>`
  type(rgb_color),parameter :: grey49               =rgb_color(125,125,125)  !:raw-html:`<span class="colored-square" style="background-color: #7D7D7D;"></span>`
  type(rgb_color),parameter :: gray50               =rgb_color(127,127,127)  !:raw-html:`<span class="colored-square" style="background-color: #7F7F7F;"></span>`
  type(rgb_color),parameter :: grey50               =rgb_color(127,127,127)  !:raw-html:`<span class="colored-square" style="background-color: #7F7F7F;"></span>`
  type(rgb_color),parameter :: gray51               =rgb_color(130,130,130)  !:raw-html:`<span class="colored-square" style="background-color: #828282;"></span>`
  type(rgb_color),parameter :: grey51               =rgb_color(130,130,130)  !:raw-html:`<span class="colored-square" style="background-color: #828282;"></span>`
  type(rgb_color),parameter :: gray52               =rgb_color(133,133,133)  !:raw-html:`<span class="colored-square" style="background-color: #858585;"></span>`
  type(rgb_color),parameter :: grey52               =rgb_color(133,133,133)  !:raw-html:`<span class="colored-square" style="background-color: #858585;"></span>`
  type(rgb_color),parameter :: gray53               =rgb_color(135,135,135)  !:raw-html:`<span class="colored-square" style="background-color: #878787;"></span>`
  type(rgb_color),parameter :: grey53               =rgb_color(135,135,135)  !:raw-html:`<span class="colored-square" style="background-color: #878787;"></span>`
  type(rgb_color),parameter :: gray54               =rgb_color(138,138,138)  !:raw-html:`<span class="colored-square" style="background-color: #8A8A8A;"></span>`
  type(rgb_color),parameter :: grey54               =rgb_color(138,138,138)  !:raw-html:`<span class="colored-square" style="background-color: #8A8A8A;"></span>`
  type(rgb_color),parameter :: gray55               =rgb_color(140,140,140)  !:raw-html:`<span class="colored-square" style="background-color: #8C8C8C;"></span>`
  type(rgb_color),parameter :: grey55               =rgb_color(140,140,140)  !:raw-html:`<span class="colored-square" style="background-color: #8C8C8C;"></span>`
  type(rgb_color),parameter :: gray56               =rgb_color(143,143,143)  !:raw-html:`<span class="colored-square" style="background-color: #8F8F8F;"></span>`
  type(rgb_color),parameter :: grey56               =rgb_color(143,143,143)  !:raw-html:`<span class="colored-square" style="background-color: #8F8F8F;"></span>`
  type(rgb_color),parameter :: gray57               =rgb_color(145,145,145)  !:raw-html:`<span class="colored-square" style="background-color: #919191;"></span>`
  type(rgb_color),parameter :: grey57               =rgb_color(145,145,145)  !:raw-html:`<span class="colored-square" style="background-color: #919191;"></span>`
  type(rgb_color),parameter :: gray58               =rgb_color(148,148,148)  !:raw-html:`<span class="colored-square" style="background-color: #949494;"></span>`
  type(rgb_color),parameter :: grey58               =rgb_color(148,148,148)  !:raw-html:`<span class="colored-square" style="background-color: #949494;"></span>`
  type(rgb_color),parameter :: gray59               =rgb_color(150,150,150)  !:raw-html:`<span class="colored-square" style="background-color: #969696;"></span>`
  type(rgb_color),parameter :: grey59               =rgb_color(150,150,150)  !:raw-html:`<span class="colored-square" style="background-color: #969696;"></span>`
  type(rgb_color),parameter :: gray60               =rgb_color(153,153,153)  !:raw-html:`<span class="colored-square" style="background-color: #999999;"></span>`
  type(rgb_color),parameter :: grey60               =rgb_color(153,153,153)  !:raw-html:`<span class="colored-square" style="background-color: #999999;"></span>`
  type(rgb_color),parameter :: gray61               =rgb_color(156,156,156)  !:raw-html:`<span class="colored-square" style="background-color: #9C9C9C;"></span>`
  type(rgb_color),parameter :: grey61               =rgb_color(156,156,156)  !:raw-html:`<span class="colored-square" style="background-color: #9C9C9C;"></span>`
  type(rgb_color),parameter :: gray62               =rgb_color(158,158,158)  !:raw-html:`<span class="colored-square" style="background-color: #9E9E9E;"></span>`
  type(rgb_color),parameter :: grey62               =rgb_color(158,158,158)  !:raw-html:`<span class="colored-square" style="background-color: #9E9E9E;"></span>`
  type(rgb_color),parameter :: gray63               =rgb_color(161,161,161)  !:raw-html:`<span class="colored-square" style="background-color: #A1A1A1;"></span>`
  type(rgb_color),parameter :: grey63               =rgb_color(161,161,161)  !:raw-html:`<span class="colored-square" style="background-color: #A1A1A1;"></span>`
  type(rgb_color),parameter :: gray64               =rgb_color(163,163,163)  !:raw-html:`<span class="colored-square" style="background-color: #A3A3A3;"></span>`
  type(rgb_color),parameter :: grey64               =rgb_color(163,163,163)  !:raw-html:`<span class="colored-square" style="background-color: #A3A3A3;"></span>`
  type(rgb_color),parameter :: gray65               =rgb_color(166,166,166)  !:raw-html:`<span class="colored-square" style="background-color: #A6A6A6;"></span>`
  type(rgb_color),parameter :: grey65               =rgb_color(166,166,166)  !:raw-html:`<span class="colored-square" style="background-color: #A6A6A6;"></span>`
  type(rgb_color),parameter :: gray66               =rgb_color(168,168,168)  !:raw-html:`<span class="colored-square" style="background-color: #A8A8A8;"></span>`
  type(rgb_color),parameter :: grey66               =rgb_color(168,168,168)  !:raw-html:`<span class="colored-square" style="background-color: #A8A8A8;"></span>`
  type(rgb_color),parameter :: gray67               =rgb_color(171,171,171)  !:raw-html:`<span class="colored-square" style="background-color: #ABABAB;"></span>`
  type(rgb_color),parameter :: grey67               =rgb_color(171,171,171)  !:raw-html:`<span class="colored-square" style="background-color: #ABABAB;"></span>`
  type(rgb_color),parameter :: gray68               =rgb_color(173,173,173)  !:raw-html:`<span class="colored-square" style="background-color: #ADADAD;"></span>`
  type(rgb_color),parameter :: grey68               =rgb_color(173,173,173)  !:raw-html:`<span class="colored-square" style="background-color: #ADADAD;"></span>`
  type(rgb_color),parameter :: gray69               =rgb_color(176,176,176)  !:raw-html:`<span class="colored-square" style="background-color: #B0B0B0;"></span>`
  type(rgb_color),parameter :: grey69               =rgb_color(176,176,176)  !:raw-html:`<span class="colored-square" style="background-color: #B0B0B0;"></span>`
  type(rgb_color),parameter :: gray70               =rgb_color(179,179,179)  !:raw-html:`<span class="colored-square" style="background-color: #B3B3B3;"></span>`
  type(rgb_color),parameter :: grey70               =rgb_color(179,179,179)  !:raw-html:`<span class="colored-square" style="background-color: #B3B3B3;"></span>`
  type(rgb_color),parameter :: gray71               =rgb_color(181,181,181)  !:raw-html:`<span class="colored-square" style="background-color: #B5B5B5;"></span>`
  type(rgb_color),parameter :: grey71               =rgb_color(181,181,181)  !:raw-html:`<span class="colored-square" style="background-color: #B5B5B5;"></span>`
  type(rgb_color),parameter :: gray72               =rgb_color(184,184,184)  !:raw-html:`<span class="colored-square" style="background-color: #B8B8B8;"></span>`
  type(rgb_color),parameter :: grey72               =rgb_color(184,184,184)  !:raw-html:`<span class="colored-square" style="background-color: #B8B8B8;"></span>`
  type(rgb_color),parameter :: gray73               =rgb_color(186,186,186)  !:raw-html:`<span class="colored-square" style="background-color: #BABABA;"></span>`
  type(rgb_color),parameter :: grey73               =rgb_color(186,186,186)  !:raw-html:`<span class="colored-square" style="background-color: #BABABA;"></span>`
  type(rgb_color),parameter :: gray74               =rgb_color(189,189,189)  !:raw-html:`<span class="colored-square" style="background-color: #BDBDBD;"></span>`
  type(rgb_color),parameter :: grey74               =rgb_color(189,189,189)  !:raw-html:`<span class="colored-square" style="background-color: #BDBDBD;"></span>`
  type(rgb_color),parameter :: gray75               =rgb_color(191,191,191)  !:raw-html:`<span class="colored-square" style="background-color: #BFBFBF;"></span>`
  type(rgb_color),parameter :: grey75               =rgb_color(191,191,191)  !:raw-html:`<span class="colored-square" style="background-color: #BFBFBF;"></span>`
  type(rgb_color),parameter :: gray76               =rgb_color(194,194,194)  !:raw-html:`<span class="colored-square" style="background-color: #C2C2C2;"></span>`
  type(rgb_color),parameter :: grey76               =rgb_color(194,194,194)  !:raw-html:`<span class="colored-square" style="background-color: #C2C2C2;"></span>`
  type(rgb_color),parameter :: gray77               =rgb_color(196,196,196)  !:raw-html:`<span class="colored-square" style="background-color: #C4C4C4;"></span>`
  type(rgb_color),parameter :: grey77               =rgb_color(196,196,196)  !:raw-html:`<span class="colored-square" style="background-color: #C4C4C4;"></span>`
  type(rgb_color),parameter :: gray78               =rgb_color(199,199,199)  !:raw-html:`<span class="colored-square" style="background-color: #C7C7C7;"></span>`
  type(rgb_color),parameter :: grey78               =rgb_color(199,199,199)  !:raw-html:`<span class="colored-square" style="background-color: #C7C7C7;"></span>`
  type(rgb_color),parameter :: gray79               =rgb_color(201,201,201)  !:raw-html:`<span class="colored-square" style="background-color: #C9C9C9;"></span>`
  type(rgb_color),parameter :: grey79               =rgb_color(201,201,201)  !:raw-html:`<span class="colored-square" style="background-color: #C9C9C9;"></span>`
  type(rgb_color),parameter :: gray80               =rgb_color(204,204,204)  !:raw-html:`<span class="colored-square" style="background-color: #CCCCCC;"></span>`
  type(rgb_color),parameter :: grey80               =rgb_color(204,204,204)  !:raw-html:`<span class="colored-square" style="background-color: #CCCCCC;"></span>`
  type(rgb_color),parameter :: gray81               =rgb_color(207,207,207)  !:raw-html:`<span class="colored-square" style="background-color: #CFCFCF;"></span>`
  type(rgb_color),parameter :: grey81               =rgb_color(207,207,207)  !:raw-html:`<span class="colored-square" style="background-color: #CFCFCF;"></span>`
  type(rgb_color),parameter :: gray82               =rgb_color(209,209,209)  !:raw-html:`<span class="colored-square" style="background-color: #D1D1D1;"></span>`
  type(rgb_color),parameter :: grey82               =rgb_color(209,209,209)  !:raw-html:`<span class="colored-square" style="background-color: #D1D1D1;"></span>`
  type(rgb_color),parameter :: gray83               =rgb_color(212,212,212)  !:raw-html:`<span class="colored-square" style="background-color: #D4D4D4;"></span>`
  type(rgb_color),parameter :: grey83               =rgb_color(212,212,212)  !:raw-html:`<span class="colored-square" style="background-color: #D4D4D4;"></span>`
  type(rgb_color),parameter :: gray84               =rgb_color(214,214,214)  !:raw-html:`<span class="colored-square" style="background-color: #D6D6D6;"></span>`
  type(rgb_color),parameter :: grey84               =rgb_color(214,214,214)  !:raw-html:`<span class="colored-square" style="background-color: #D6D6D6;"></span>`
  type(rgb_color),parameter :: gray85               =rgb_color(217,217,217)  !:raw-html:`<span class="colored-square" style="background-color: #D9D9D9;"></span>`
  type(rgb_color),parameter :: grey85               =rgb_color(217,217,217)  !:raw-html:`<span class="colored-square" style="background-color: #D9D9D9;"></span>`
  type(rgb_color),parameter :: gray86               =rgb_color(219,219,219)  !:raw-html:`<span class="colored-square" style="background-color: #DBDBDB;"></span>`
  type(rgb_color),parameter :: grey86               =rgb_color(219,219,219)  !:raw-html:`<span class="colored-square" style="background-color: #DBDBDB;"></span>`
  type(rgb_color),parameter :: gray87               =rgb_color(222,222,222)  !:raw-html:`<span class="colored-square" style="background-color: #DEDEDE;"></span>`
  type(rgb_color),parameter :: grey87               =rgb_color(222,222,222)  !:raw-html:`<span class="colored-square" style="background-color: #DEDEDE;"></span>`
  type(rgb_color),parameter :: gray88               =rgb_color(224,224,224)  !:raw-html:`<span class="colored-square" style="background-color: #E0E0E0;"></span>`
  type(rgb_color),parameter :: grey88               =rgb_color(224,224,224)  !:raw-html:`<span class="colored-square" style="background-color: #E0E0E0;"></span>`
  type(rgb_color),parameter :: gray89               =rgb_color(227,227,227)  !:raw-html:`<span class="colored-square" style="background-color: #E3E3E3;"></span>`
  type(rgb_color),parameter :: grey89               =rgb_color(227,227,227)  !:raw-html:`<span class="colored-square" style="background-color: #E3E3E3;"></span>`
  type(rgb_color),parameter :: gray90               =rgb_color(229,229,229)  !:raw-html:`<span class="colored-square" style="background-color: #E5E5E5;"></span>`
  type(rgb_color),parameter :: grey90               =rgb_color(229,229,229)  !:raw-html:`<span class="colored-square" style="background-color: #E5E5E5;"></span>`
  type(rgb_color),parameter :: gray91               =rgb_color(232,232,232)  !:raw-html:`<span class="colored-square" style="background-color: #E8E8E8;"></span>`
  type(rgb_color),parameter :: grey91               =rgb_color(232,232,232)  !:raw-html:`<span class="colored-square" style="background-color: #E8E8E8;"></span>`
  type(rgb_color),parameter :: gray92               =rgb_color(235,235,235)  !:raw-html:`<span class="colored-square" style="background-color: #EBEBEB;"></span>`
  type(rgb_color),parameter :: grey92               =rgb_color(235,235,235)  !:raw-html:`<span class="colored-square" style="background-color: #EBEBEB;"></span>`
  type(rgb_color),parameter :: gray93               =rgb_color(237,237,237)  !:raw-html:`<span class="colored-square" style="background-color: #EDEDED;"></span>`
  type(rgb_color),parameter :: grey93               =rgb_color(237,237,237)  !:raw-html:`<span class="colored-square" style="background-color: #EDEDED;"></span>`
  type(rgb_color),parameter :: gray94               =rgb_color(240,240,240)  !:raw-html:`<span class="colored-square" style="background-color: #F0F0F0;"></span>`
  type(rgb_color),parameter :: grey94               =rgb_color(240,240,240)  !:raw-html:`<span class="colored-square" style="background-color: #F0F0F0;"></span>`
  type(rgb_color),parameter :: gray95               =rgb_color(242,242,242)  !:raw-html:`<span class="colored-square" style="background-color: #F2F2F2;"></span>`
  type(rgb_color),parameter :: grey95               =rgb_color(242,242,242)  !:raw-html:`<span class="colored-square" style="background-color: #F2F2F2;"></span>`
  type(rgb_color),parameter :: gray96               =rgb_color(245,245,245)  !:raw-html:`<span class="colored-square" style="background-color: #F5F5F5;"></span>`
  type(rgb_color),parameter :: grey96               =rgb_color(245,245,245)  !:raw-html:`<span class="colored-square" style="background-color: #F5F5F5;"></span>`
  type(rgb_color),parameter :: gray97               =rgb_color(247,247,247)  !:raw-html:`<span class="colored-square" style="background-color: #F7F7F7;"></span>`
  type(rgb_color),parameter :: grey97               =rgb_color(247,247,247)  !:raw-html:`<span class="colored-square" style="background-color: #F7F7F7;"></span>`
  type(rgb_color),parameter :: gray98               =rgb_color(250,250,250)  !:raw-html:`<span class="colored-square" style="background-color: #FAFAFA;"></span>`
  type(rgb_color),parameter :: grey98               =rgb_color(250,250,250)  !:raw-html:`<span class="colored-square" style="background-color: #FAFAFA;"></span>`
  type(rgb_color),parameter :: gray99               =rgb_color(252,252,252)  !:raw-html:`<span class="colored-square" style="background-color: #FCFCFC;"></span>`
  type(rgb_color),parameter :: grey99               =rgb_color(252,252,252)  !:raw-html:`<span class="colored-square" style="background-color: #FCFCFC;"></span>`
  type(rgb_color),parameter :: gray100              =rgb_color(255,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFFF;"></span>`
  type(rgb_color),parameter :: grey100              =rgb_color(255,255,255)  !:raw-html:`<span class="colored-square" style="background-color: #FFFFFF;"></span>`
  type(rgb_color),parameter :: darkgrey             =rgb_color(169,169,169)  !:raw-html:`<span class="colored-square" style="background-color: #A9A9A9;"></span>`
  type(rgb_color),parameter :: darkgray             =rgb_color(169,169,169)  !:raw-html:`<span class="colored-square" style="background-color: #A9A9A9;"></span>`
  type(rgb_color),parameter :: darkblue             =rgb_color(0,0,139)  !:raw-html:`<span class="colored-square" style="background-color: #00008B;"></span>`
  type(rgb_color),parameter :: darkcyan             =rgb_color(0,139,139)  !:raw-html:`<span class="colored-square" style="background-color: #008B8B;"></span>`
  type(rgb_color),parameter :: darkmagenta          =rgb_color(139,0,139)  !:raw-html:`<span class="colored-square" style="background-color: #8B008B;"></span>`
  type(rgb_color),parameter :: darkred              =rgb_color(139,0,0)  !:raw-html:`<span class="colored-square" style="background-color: #8B0000;"></span>`
  type(rgb_color),parameter :: lightgreen           =rgb_color(144,238,144)  !:raw-html:`<span class="colored-square" style="background-color: #90EE90;"></span>`


  interface assignment(=)
  !
  !Overloaded assignment operator for :f:type:`rgb_color` 
  !types using :f:func:`equal_colors`
  !
     module procedure equal_colors
  end interface assignment(=)

  interface operator(+)
  !
  !Overloaded addition operator for :f:type:`rgb_color` 
  !types using :f:func:`add_colors`
  !
     module procedure add_colors
  end interface operator(+)

  interface operator(-)
  !
  !Overloaded subtraction operator for :f:type:`rgb_color` 
  !types using :f:func:`subtract_colors`
  !
     module procedure subtract_colors
  end interface operator(-)

  interface operator(*)
  !
  !Overloaded multiplication operator for :f:type:`rgb_color` 
  !types using :f:func:`scalar_left_colors`
  !
     module procedure scalar_left_color
  end interface operator(*)

  interface operator(.dot.)
  !
  !Overloaded dot operator for :f:type:`rgb_color` 
  !types using :f:func:`dot_scalar_colors`
  !
     module procedure dot_scalar_colors
  end interface operator(.dot.)

contains

  function rgb(c) result(num)
  !
  !Given a :f:type:`rgb_color` :f:var:`c` :code:`=[r,g,b]`, returns 
  !an univocally determined integer number :code:`65536*r+256*g+b`
  !
    type(rgb_color),intent(in) :: c !the input :f:type:`rgb_color`
    integer :: num                  !the output number
    num = int(c%r)*65536 + int(c%g)*256 + int(c%b)
  end function rgb

  elemental subroutine equal_colors(C1,C2)
  !
  !Given a :f:type:`rgb_color` :f:var:`c2` :code:`=[r,g,b]`, assigns 
  !its value to a a :f:type:`rgb_color` :f:var:`c1`
  !
    type(rgb_color),intent(in)    :: C2 !:f:type:`rgb_color` on the right-hand side of the equality
    type(rgb_color),intent(inout) :: C1 !:f:type:`rgb_color` on the left-hand side of the equality
    C1%r = C2%r
    C1%g = C2%g
    C1%b = C2%b
  end subroutine equal_colors

  elemental function add_colors(c1,c2) result(c)
  !
  !Given two :f:type:`rgb_color` :code:`[r1,g1,b1]` and :code:`[r2,g2,b2]`
  !returns :code:`[r1+r2,g1+g2,b1+b2]`.
  !
    type(rgb_color),intent(in) :: c1  !First color to sum
    type(rgb_color),intent(in) :: c2  !Second color to sum
    type(rgb_color)            :: c   !Resulting color
    c%r = c1%r + c2%r
    c%g = c1%g + c2%g
    c%b = c1%b + c2%b
  end function add_colors

  elemental function subtract_colors(c1,c2) result(c)
  !
  !Given two :f:type:`rgb_color` :code:`[r1,g1,b1]` and :code:`[r2,g2,b2]`
  !returns :code:`[r1-r2,g1-g2,b1-b2]`.
  !
    type(rgb_color),intent(in) :: c1  !Minuend
    type(rgb_color),intent(in) :: c2  !Subtrahent
    type(rgb_color)            :: c   !Result
    c%r = c1%r - c2%r
    c%g = c1%g - c2%g
    c%b = c1%b - c2%b
  end function subtract_colors

  elemental function scalar_left_color(k,cin) result(cout)
  !
  !Given a real :f:var:`k` and a :f:type:`rgb_color` :code:`[r,g,b]`
  !returns :f:type:`rgb_color` :code:`[k*r,k*g,k*b]`.
  !
    real(8),intent(in)         :: k    !scalar term  (multiplier)
    type(rgb_color),intent(in) :: cin  !rgb color (multiplicand)
    type(rgb_color)            :: cout !result
    cout%r = k*cin%r
    cout%g = k*cin%g
    cout%b = k*cin%b
  end function scalar_left_color

  elemental function scalar_right_color(k,cin) result(cout)
  !
  !Given an :f:type:`rgb_color` :code:`[r,g,b]` and a real :f:var:`k`
  !returns :f:type:`rgb_color` :code:`[k*r,k*g,k*b]`.
  !
    real(8),intent(in)         :: k    !scalar term  (multiplier)
    type(rgb_color),intent(in) :: cin  !rgb color (multiplicand)
    type(rgb_color)            :: cout !result
    cout%r = cin%r*k
    cout%g = cin%g*k
    cout%b = cin%b*k
  end function scalar_right_color

  function dot_scalar_colors(v,cin) result(cout)
  !
  !Given an array of real :f:var:`v` and an array of :f:type:`rgb_color` :f:var:`cin`
  !returns a :f:type:`rgb_color` :f:var:`cout` = :math:`\sum_i vin(i)*cin(i)` 
  !
    real(8),dimension(:),intent(in)               :: v    !array of reals
    type(rgb_color),dimension(size(v)),intent(in) :: cin  !array of :f:type:`rgb_color`
    type(rgb_color)                               :: cout !resulting :f:type:`rgb_color`
    integer :: i
    cout=rgb_color(0,0,0)
    do i=1,size(v)
       cout%r = cout%r + v(i)*cin(i)%r
       cout%g = cout%g + v(i)*cin(i)%g
       cout%b = cout%b + v(i)*cin(i)%b
    enddo
  end function dot_scalar_colors

  function pick_color(string) result(crgb)
  !
  !Given an string, picks a :f:type:`rgb_color` among the default ones 
  !:f:var:`black`, :f:var:`blue`, :f:var:`cyan`, :f:var:`green`, :f:var:`magenta`,
  !:f:var:`orange`, :f:var:`red`, :f:var:`white`, :f:var:`yellow`
  !
    character(len=*)                                     :: string  ! name of the color
    type(rgb_color)                                      :: crgb    ! resulting :f:var:`rgb_color`
    character(len=len_trim(trim(adjustl(trim(string))))) :: color_name
    color_name=trim(adjustl(trim(string)))
    select case(color_name)
    case("black")
       crgb=black
    case("red")
       crgb=red
    case("green")
       crgb=green
    case("orange")
       crgb=orange
    case("blue")
       crgb=blue
    case("yellow")
       crgb=yellow
    case("cyan")
       crgb=cyan
    case("magenta")
       crgb=magenta
    case("white")
       crgb=white
    case default
       print*,"pick_color: color name ",color_name," does not exist"
       print*,"set color to black"
       crgb=black
    end select
  end function pick_color

END MODULE SF_COLORS
