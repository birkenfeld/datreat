
MODULE gr_interface
! compile: gfortran gr_interface.f90 ....f90   -Wl,-rpath,/usr/local/gr/lib  -L/usr/local/gr/lib -lGR
 
implicit none



  interface gr_fmt
     module procedure  gr_fmt_int, gr_fmt_dbl
  end interface

  character(len=80) :: PDF_VIEWER_COMMAND = "open -a preview "   ! mac
  character(len=80) :: LAST_DTR_PLOT


  character(len=4) :: czero
  integer          :: izero = 0
  Equivalence( izero, czero )

! Parameters for GR-framework, colors are still missing
  integer, parameter :: WS_TYPE_XTERM              = 211   !  xterm
  integer, parameter :: WS_TYPE_EPS                = 62    !  eps
  integer, parameter :: WS_TYPE_PDF                = 101   !  pdf
  integer, parameter :: WS_TYPE_PDF_COMPRESSED     = 102   !  pdf

  double precision, parameter   :: DEVICE_LOWER_LEFT_X          = 0.00d0
  double precision, parameter   :: DEVICE_LOWER_LEFT_Y          = 0.00d0
  double precision, parameter   :: DEVICE_UPPER_RIGHT_X         = 0.32d0
  double precision, parameter   :: DEVICE_UPPER_RIGHT_Y         =  DEVICE_UPPER_RIGHT_X / 1.618d0

  double precision, parameter   :: DEFAULT_WC_MARGIN             = 0.1d0

  double precision              :: AXIS_LINEWIDTH                = 1d0
  double precision              :: XLEG_DISTANCE                 = 0.075d0
  double precision              :: YLEG_DISTANCE                 = 0.14d0
  double precision              :: TLEG_DISTANCE_X               = 0.0d0
  double precision              :: TLEG_DISTANCE_Y               = 0.08d0

                         
  integer, parameter :: MARKERTYPE_DOT             =  1    !    Smallest displayable dot
  integer, parameter :: MARKERTYPE_PLUS            =  2    !    Plus sign
  integer, parameter :: MARKERTYPE_ASTERISK        =  3    !    Asterisk
  integer, parameter :: MARKERTYPE_CIRCLE          =  4    !    Hollow circle
  integer, parameter :: MARKERTYPE_DIAGONAL_CROSS  =  5    !    Diagonal cross
  integer, parameter :: MARKERTYPE_SOLID_CIRCLE    = -1    !    Filled circle
  integer, parameter :: MARKERTYPE_TRIANGLE_UP     = -2    !    Hollow triangle pointing upward
  integer, parameter :: MARKERTYPE_SOLID_TRI_UP    = -3    !    Filled triangle pointing upward
  integer, parameter :: MARKERTYPE_TRIANGLE_DOWN   = -4    !    Hollow triangle pointing downward
  integer, parameter :: MARKERTYPE_SOLID_TRI_DOWN  = -5    !    Filled triangle pointing downward
  integer, parameter :: MARKERTYPE_SQUARE          = -6    !    Hollow square
  integer, parameter :: MARKERTYPE_SOLID_SQUARE    = -7    !    Filled square
  integer, parameter :: MARKERTYPE_BOWTIE          = -8    !    Hollow bowtie
  integer, parameter :: MARKERTYPE_SOLID_BOWTIE    = -9    !    Filled bowtie
  integer, parameter :: MARKERTYPE_HGLASS          = -10   !    Hollow hourglass
  integer, parameter :: MARKERTYPE_SOLID_HGLASS    = -11   !    Filled hourglass
  integer, parameter :: MARKERTYPE_DIAMOND         = -12   !    Hollow diamond
  integer, parameter :: MARKERTYPE_SOLID_DIAMOND   = -13   !    Filled Diamond
  integer, parameter :: MARKERTYPE_STAR            = -14   !    Hollow star
  integer, parameter :: MARKERTYPE_SOLID_STAR      = -15   !    Filled Star
  integer, parameter :: MARKERTYPE_TRI_UP_DOWN     = -16   !    Hollow triangles pointing up and down overlaid
  integer, parameter :: MARKERTYPE_SOLID_TRI_RIGHT = -17   !    Filled triangle point right
  integer, parameter :: MARKERTYPE_SOLID_TRI_LEFT  = -18   !    Filled triangle pointing left
  integer, parameter :: MARKERTYPE_HOLLOW_PLUS     = -19   !    Hollow plus sign
  integer, parameter :: MARKERTYPE_SOLID_PLUS      = -20   !    Solid plus sign
  integer, parameter :: MARKERTYPE_PENTAGON        = -21   !    Pentagon
  integer, parameter :: MARKERTYPE_HEXAGON         = -22   !    Hexagon
  integer, parameter :: MARKERTYPE_HEPTAGON        = -23   !    Heptagon
  integer, parameter :: MARKERTYPE_OCTAGON         = -24   !    Octagon
  integer, parameter :: MARKERTYPE_STAR_4          = -25   !    4-pointed star
  integer, parameter :: MARKERTYPE_STAR_5          = -26   !    5-pointed star (pentagram)
  integer, parameter :: MARKERTYPE_STAR_6          = -27   !    6-pointed star (hexagram)
  integer, parameter :: MARKERTYPE_STAR_7          = -28   !    7-pointed star (heptagram)
  integer, parameter :: MARKERTYPE_STAR_8          = -29   !    8-pointed star (octagram)
  integer, parameter :: MARKERTYPE_VLINE           = -30   !    verical line
  integer, parameter :: MARKERTYPE_HLINE           = -31   !    horizontal line
  integer, parameter :: MARKERTYPE_OMARK           = -32   !    o-mark

  integer, parameter :: MARKERTYPE(37) = [                        &
                                       MARKERTYPE_DOT            ,&  
                                       MARKERTYPE_CIRCLE         ,& 
                                       MARKERTYPE_SQUARE         ,&
                                       MARKERTYPE_DIAMOND        ,&
                                       MARKERTYPE_TRIANGLE_UP    ,&
                                       MARKERTYPE_TRIANGLE_DOWN  ,&
                                       MARKERTYPE_BOWTIE         ,&
                                       MARKERTYPE_SOLID_TRI_UP   ,&
                                       MARKERTYPE_SOLID_TRI_DOWN ,&
                                       MARKERTYPE_SOLID_SQUARE   ,&
                                       MARKERTYPE_PLUS           ,&
                                       MARKERTYPE_ASTERISK       ,&
                                       MARKERTYPE_DIAGONAL_CROSS ,&
                                       MARKERTYPE_SOLID_CIRCLE   ,&
                                       MARKERTYPE_SOLID_BOWTIE   ,&
                                       MARKERTYPE_HGLASS         ,&
                                       MARKERTYPE_SOLID_HGLASS   ,&
                                       MARKERTYPE_SOLID_DIAMOND  ,&
                                       MARKERTYPE_STAR           ,&
                                       MARKERTYPE_SOLID_STAR     ,&
                                       MARKERTYPE_TRI_UP_DOWN    ,&
                                       MARKERTYPE_SOLID_TRI_RIGHT,&
                                       MARKERTYPE_SOLID_TRI_LEFT ,&
                                       MARKERTYPE_HOLLOW_PLUS    ,&
                                       MARKERTYPE_SOLID_PLUS     ,&
                                       MARKERTYPE_PENTAGON       ,&
                                       MARKERTYPE_HEXAGON        ,&
                                       MARKERTYPE_HEPTAGON       ,&
                                       MARKERTYPE_OCTAGON        ,&
                                       MARKERTYPE_STAR_4         ,&
                                       MARKERTYPE_STAR_5         ,&
                                       MARKERTYPE_STAR_6         ,&
                                       MARKERTYPE_STAR_7         ,&
                                       MARKERTYPE_STAR_8         ,&
                                       MARKERTYPE_VLINE          ,&
                                       MARKERTYPE_HLINE          ,&
                                       MARKERTYPE_OMARK           &
                                         ]


  integer, parameter :: LINETYPE_SOLID           =  1  !   Solid line
  integer, parameter :: LINETYPE_DASHED          =  2  !   Dashed line
  integer, parameter :: LINETYPE_DOTTED          =  3  !   Dotted line
  integer, parameter :: LINETYPE_DASHED_DOTTED   =  4  !   Dashed-dotted line
  integer, parameter :: LINETYPE_DASH_2_DOT      = -1  !   Sequence of one dash followed by two dots
  integer, parameter :: LINETYPE_DASH_3_DOT      = -2  !   Sequence of one dash followed by three dots
  integer, parameter :: LINETYPE_LONG_DASH       = -3  !   Sequence of long dashes
  integer, parameter :: LINETYPE_LONG_SHORT_DASH = -4  !   Sequence of a long dash followed by a short dash
  integer, parameter :: LINETYPE_SPACED_DASH     = -5  !   Sequence of dashes double spaced
  integer, parameter :: LINETYPE_SPACED_DOT      = -6  !   Sequence of dots double spaced
  integer, parameter :: LINETYPE_DOUBLE_DOT      = -7  !   Sequence of pairs of dots
  integer, parameter :: LINETYPE_TRIPLE_DOT      = -8  !   Sequence of groups of three dots

  
  integer, parameter :: OPTION_LINEAR  = 0 ! ?? 
  integer, parameter :: OPTION_X_LOG   = 1 ! ?? 
  integer, parameter :: OPTION_Y_LOG   = 2
  integer, parameter :: OPTION_XY_LOG  = 3
  integer, parameter :: OPTION_FLIP_X  = 4
  integer, parameter :: OPTION_FLIP_Y  = 5
  integer, parameter :: OPTION_FLIP_Z  = 6 ! ?? 


! Define the current direction in which subsequent text will be drawn.

  integer, parameter  :: TEXT_PATH_RIGHT  =  0   !  left-to-right
  integer, parameter  :: TEXT_PATH_LEFT   =  1   !  right-to-left
  integer, parameter  :: TEXT_PATH_UP     =  2   !  downside-up
  integer, parameter  :: TEXT_PATH_DOWN   =  3   !  upside-down



! Set the current horizontal and vertical alignment for text.

   integer, parameter  :: TEXT_HALIGN_NORMAL = 0  !    
   integer, parameter  :: TEXT_HALIGN_LEFT   = 1  !   Left justify
   integer, parameter  :: TEXT_HALIGN_CENTER = 2  !   Center justify
   integer, parameter  :: TEXT_HALIGN_RIGHT  = 3  !   Right justify
   integer, parameter  :: TEXT_VALIGN_NORMAL = 0  !    
   integer, parameter  :: TEXT_VALIGN_TOP    = 1  !   Align with the top of the characters
   integer, parameter  :: TEXT_VALIGN_CAP    = 2  !   Aligned with the cap of the characters
   integer, parameter  :: TEXT_VALIGN_HALF   = 3  !   Aligned with the half line of the characters
   integer, parameter  :: TEXT_VALIGN_BASE   = 4  !   Aligned with the base line of the characters
   integer, parameter  :: TEXT_VALIGN_BOTTOM = 5  !   Aligned with the bottom line of the characters

! text fonts
   integer, parameter  :: FONT_TIMES_ROMAN                    = 101
   integer, parameter  :: FONT_TIMES_ITALIC                   = 102
   integer, parameter  :: FONT_TIMES_BOLD                     = 103
   integer, parameter  :: FONT_TIMES_BOLDITALIC               = 104
   integer, parameter  :: FONT_HELVETICA                      = 105
   integer, parameter  :: FONT_HELVETICA_OBLIQUE              = 106
   integer, parameter  :: FONT_HELVETICA_BOLD                 = 107
   integer, parameter  :: FONT_HELVETICA_BOLDOBLIQUE          = 108
   integer, parameter  :: FONT_COURIER                        = 109
   integer, parameter  :: FONT_COURIER_OBLIQUE                = 110
   integer, parameter  :: FONT_COURIER_BOLD                   = 111
   integer, parameter  :: FONT_COURIER_BOLDOBLIQUE            = 112
   integer, parameter  :: FONT_SYMBOL                         = 113
   integer, parameter  :: FONT_BOOKMAN_LIGHT                  = 114
   integer, parameter  :: FONT_BOOKMAN_LIGHTITALIC            = 115
   integer, parameter  :: FONT_BOOKMAN_DEMI                   = 116
   integer, parameter  :: FONT_BOOKMAN_DEMIITALIC             = 117
   integer, parameter  :: FONT_NEWCENTURYSCHLBK_ROMAN         = 118
   integer, parameter  :: FONT_NEWCENTURYSCHLBK_ITALIC        = 119
   integer, parameter  :: FONT_NEWCENTURYSCHLBK_BOLD          = 120
   integer, parameter  :: FONT_NEWCENTURYSCHLBK_BOLDITALIC    = 121
   integer, parameter  :: FONT_AVANTGARDE_BOOK                = 122
   integer, parameter  :: FONT_AVANTGARDE_BOOKOBLIQUE         = 123
   integer, parameter  :: FONT_AVANTGARDE_DEMI                = 124
   integer, parameter  :: FONT_AVANTGARDE_DEMIOBLIQUE         = 125
   integer, parameter  :: FONT_PALATINO_ROMAN                 = 126
   integer, parameter  :: FONT_PALATINO_ITALIC                = 127
   integer, parameter  :: FONT_PALATINO_BOLD                  = 128
   integer, parameter  :: FONT_PALATINO_BOLDITALIC            = 129
   integer, parameter  :: FONT_ZAPFCHANCERY_MEDIUMITALIC      = 130
   integer, parameter  :: FONT_ZAPFDINGBATS                   = 131

! The available text precisions are:
   
   integer, parameter  :: TEXT_PRECISION_STRING  = 0  !  String precision (higher quality)
   integer, parameter  :: TEXT_PRECISION_CHAR    = 1  ! Character precision (medium quality)
   integer, parameter  :: TEXT_PRECISION_STROKE  = 2  !  Stroke precision (lower quality)


! void gr_setfillintstyle(int style)
! Set the fill area interior style to be used for fill areas.

   integer, parameter  :: HOLLOW   = 0  !   No filling. Just draw the bounding polyline
   integer, parameter  :: SOLID    = 1  !   Fill the interior of the polygon using the fill color index
   integer, parameter  :: PATTERN  = 2  !   Fill the interior of the polygon using the style index as a pattern index
   integer, parameter  :: HATCH    = 3  !   Fill the interior of the polygon using the style index as a cross-hatched style

! COLORS
   integer, parameter  :: GR_WHITE   =  0
   integer, parameter  :: GR_BLACK   =  1
   integer, parameter  :: GR_RED     =  2
   integer, parameter  :: GR_GREEN   =  3
   integer, parameter  :: GR_BLUE    =  4
   integer, parameter  :: GR_CYAN    =  5
   integer, parameter  :: GR_YELLOW  =  6
   integer, parameter  :: GR_MAGENTA =  7
   integer, parameter  :: GR_GRAY      =  8
   integer, parameter  :: GR_LIGHTGRAY =  9
   integer, parameter  :: GR_DARKGRAY  = 10
 

   double precision, dimension(3), parameter  :: RGB_GRAY        = [ 0.50d0 , 0.50d0, 0.50d0 ]
   double precision, dimension(3), parameter  :: RGB_LIGHTGRAY   = [ 0.25d0 , 0.25d0, 0.25d0 ]
   double precision, dimension(3), parameter  :: RGB_DARKGRAY    = [ 0.75d0 , 0.75d0, 0.75d0 ]


! names <==> numbers from :  https://www.rapidtables.com/web/color/RGB_Color.html
   double  precision,  dimension(3),  parameter  ::  RGB_MAROON                   =  [128,0,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_RED                 =  [139,0,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BROWN                    =  [165,42,42]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_FIREBRICK                =  [178,34,34]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CRIMSON                  =  [220,20,60]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_RED                      =  [255,0,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_TOMATO                   =  [255,99,71]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CORAL                    =  [255,127,80]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_INDIAN_RED               =  [205,92,92]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_CORAL              =  [240,128,128]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_SALMON              =  [233,150,122]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SALMON                   =  [250,128,114]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_SALMON             =  [255,160,122]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ORANGE_RED               =  [255,69,0]     /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_ORANGE              =  [255,140,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ORANGE                   =  [255,165,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GOLD                     =  [255,215,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_GOLDEN_ROD          =  [184,134,11]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GOLDEN_ROD               =  [218,165,32]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PALE_GOLDEN_ROD          =  [238,232,170]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_KHAKI               =  [189,183,107]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_KHAKI                    =  [240,230,140]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_OLIVE                    =  [128,128,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_YELLOW                   =  [255,255,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_YELLOW_GREEN             =  [154,205,50]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_OLIVE_GREEN         =  [85,107,47]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_OLIVE_DRAB               =  [107,142,35]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LAWN_GREEN               =  [124,252,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CHARTREUSE               =  [127,255,0]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GREEN_YELLOW             =  [173,255,47]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_GREEN               =  [0,100,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GREEN                    =  [0,128,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_FOREST_GREEN             =  [34,139,34]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIME                     =  [0,255,0]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIME_GREEN               =  [50,205,50]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_GREEN              =  [144,238,144]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PALE_GREEN               =  [152,251,152]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_SEA_GREEN           =  [143,188,143]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_SPRING_GREEN      =  [0,250,154]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SPRING_GREEN             =  [0,255,127]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SEA_GREEN                =  [46,139,87]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_AQUA_MARINE       =  [102,205,170]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_SEA_GREEN         =  [60,179,113]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_SEA_GREEN          =  [32,178,170]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_SLATE_GRAY          =  [47,79,79]     /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_TEAL                     =  [0,128,128]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_CYAN                =  [0,139,139]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_AQUA                     =  [0,255,255]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CYAN                     =  [0,255,255]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_CYAN               =  [224,255,255]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_TURQUOISE           =  [0,206,209]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_TURQUOISE                =  [64,224,208]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_TURQUOISE         =  [72,209,204]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PALE_TURQUOISE           =  [175,238,238]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_AQUA_MARINE              =  [127,255,212]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_POWDER_BLUE              =  [176,224,230]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CADET_BLUE               =  [95,158,160]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_STEEL_BLUE               =  [70,130,180]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CORN_FLOWER_BLUE         =  [100,149,237]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DEEP_SKY_BLUE            =  [0,191,255]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DODGER_BLUE              =  [30,144,255]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_BLUE               =  [173,216,230]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SKY_BLUE                 =  [135,206,235]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_SKY_BLUE           =  [135,206,250]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MIDNIGHT_BLUE            =  [25,25,112]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_NAVY                     =  [0,0,128]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_BLUE                =  [0,0,139]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_BLUE              =  [0,0,205]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BLUE                     =  [0,0,255]      /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ROYAL_BLUE               =  [65,105,225]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BLUE_VIOLET              =  [138,43,226]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_INDIGO                   =  [75,0,130]     /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_SLATE_BLUE          =  [72,61,139]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SLATE_BLUE               =  [106,90,205]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_SLATE_BLUE        =  [123,104,238]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_PURPLE            =  [147,112,219]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_MAGENTA             =  [139,0,139]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_VIOLET              =  [148,0,211]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DARK_ORCHID              =  [153,50,204]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_ORCHID            =  [186,85,211]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PURPLE                   =  [128,0,128]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_THISTLE                  =  [216,191,216]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PLUM                     =  [221,160,221]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_VIOLET                   =  [238,130,238]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MAGENTA                  =  [255,0,255]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ORCHID                   =  [218,112,214]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MEDIUM_VIOLET_RED        =  [199,21,133]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PALE_VIOLET_RED          =  [219,112,147]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_DEEP_PINK                =  [255,20,147]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_HOT_PINK                 =  [255,105,180]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_PINK               =  [255,182,193]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PINK                     =  [255,192,203]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ANTIQUE_WHITE            =  [250,235,215]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BEIGE                    =  [245,245,220]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BISQUE                   =  [255,228,196]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BLANCHED_ALMOND          =  [255,235,205]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_WHEAT                    =  [245,222,179]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CORN_SILK                =  [255,248,220]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LEMON_CHIFFON            =  [255,250,205]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_GOLDEN_ROD_YELLOW  =  [250,250,210]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_YELLOW             =  [255,255,224]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SADDLE_BROWN             =  [139,69,19]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SIENNA                   =  [160,82,45]    /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_CHOCOLATE                =  [210,105,30]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PERU                     =  [205,133,63]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SANDY_BROWN              =  [244,164,96]   /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BURLY_WOOD               =  [222,184,135]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_TAN                      =  [210,180,140]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ROSY_BROWN               =  [188,143,143]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MOCCASIN                 =  [255,228,181]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_NAVAJO_WHITE             =  [255,222,173]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PEACH_PUFF               =  [255,218,185]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MISTY_ROSE               =  [255,228,225]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LAVENDER_BLUSH           =  [255,240,245]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LINEN                    =  [250,240,230]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_OLD_LACE                 =  [253,245,230]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_PAPAYA_WHIP              =  [255,239,213]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SEA_SHELL                =  [255,245,238]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_MINT_CREAM               =  [245,255,250]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SLATE_GRAY               =  [112,128,144]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_SLATE_GRAY         =  [119,136,153]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LIGHT_STEEL_BLUE         =  [176,196,222]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_LAVENDER                 =  [230,230,250]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_FLORAL_WHITE             =  [255,250,240]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_ALICE_BLUE               =  [240,248,255]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GHOST_WHITE              =  [248,248,255]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_HONEYDEW                 =  [240,255,240]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_IVORY                    =  [255,255,240]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_AZURE                    =  [240,255,255]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SNOW                     =  [255,250,250]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_BLACK                    =  [0,0,0]        /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_SILVER                   =  [192,192,192]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_GAINSBORO                =  [220,220,220]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_WHITE_SMOKE              =  [245,245,245]  /  255d0
   double  precision,  dimension(3),  parameter  ::  RGB_WHITE                    =  [255,255,255]  /  255d0
   

   integer, parameter           :: GRSTRLEN = 128
   integer                      :: gr_ws_id = 1   ! Workstation handle 
   character(len=GRSTRLEN)      :: gr_plotfile

   integer, parameter           :: PRIVATE_COLOR_INDEX_START = 20
   integer, private             :: next_color_index =  PRIVATE_COLOR_INDEX_START 
   integer, parameter           :: MAX_COLOR_INDEX  = 1000

   logical                      :: grinit = .false.

   double precision :: graspekt


CONTAINS                                                                       

!== gr_start =================================================================================================  
  subroutine grstart( plotfile, boundary )
     implicit none
     character(len=*), intent(in), optional  :: plotfile
     double precision, intent(in), optional  :: boundary  ! extra margin around world coordinates 


     double precision :: device_xmin = DEVICE_LOWER_LEFT_X  
     double precision :: device_xmax = DEVICE_UPPER_RIGHT_X 
     double precision :: device_ymin = DEVICE_LOWER_LEFT_Y 
     double precision :: device_ymax = DEVICE_UPPER_RIGHT_Y

     double precision :: xyboundary

     if(.not. grinit) then
       call gr_opengks()
       grinit = .true.
     endif
     if(present(plotfile) .and. len_trim(plotfile) > 0) then
        gr_plotfile = plotfile
        call gr_openws(gr_ws_id,trim(gr_plotfile)//czero,WS_TYPE_PDF_COMPRESSED)
     else
        gr_plotfile = " "
        call gr_openws(gr_ws_id,"dummy",WS_TYPE_XTERM)     ! 0=std term
     endif
     
     call gr_activatews(gr_ws_id)
                                                         
     call gr_setwsviewport(device_xmin, device_xmax, device_ymin, device_ymax) 

!! > neu
     graspekt =  (device_ymax-device_ymin)/(device_xmax-device_xmin) 
     call gr_setwswindow( 0d0, 1d0, 0d0, graspekt )
!! <
     
     if(present(boundary)) then   
       xyboundary = boundary 
     else
       xyboundary = DEFAULT_WC_MARGIN   
     endif
     call gr_setviewport(xyboundary,  & 
             (1d0-xyboundary)*graspekt, xyboundary, (1d0-xyboundary)*graspekt) 


     call gr_setmarkersize(1.0D0) 
     call gr_setmarkertype(MARKERTYPE_SOLID_CIRCLE) 
     call gr_setcharheight(0.024D0) 
     call gr_settextalign(TEXT_HALIGN_CENTER , TEXT_VALIGN_NORMAL ) 
     call gr_settextfontprec(3, 0) 


     next_color_index = PRIVATE_COLOR_INDEX_START
     call gr_setcolorrep(GR_GRAY       ,  RGB_GRAY(1)     ,  RGB_GRAY(2)     ,  RGB_GRAY(3))  
     call gr_setcolorrep(GR_LIGHTGRAY  ,  RGB_LIGHTGRAY(1),  RGB_LIGHTGRAY(2),  RGB_LIGHTGRAY(3))  
     call gr_setcolorrep(GR_DARKGRAY   ,  RGB_DARKGRAY(1) ,  RGB_DARKGRAY(2) ,  RGB_DARKGRAY(3))  

  end subroutine grstart 

!== gr_ende ================================================================================================= 
! tbd return filename (function?) depending on ws_id 
  subroutine grende(action)
     implicit none
      character(len=*), intent(in), optional :: action

      call gr_updatews() 

      if(len_trim(gr_plotfile) == 0) then
        write(*,*)"ENTER to proceed!"
        read(*,'(a)')
      endif

      call gr_deactivatews(gr_ws_id) 
      call gr_closews(gr_ws_id)   


      if(present(action)) then  
         if( .not. ( index(action,"$plot")>0 .and. len_trim(gr_plotfile)==0 ) ) then  
! do we need a "sync" in addition here ? Or a delay ?  
write(*,*)"Tgr execute:", trim(gr_string_replace(action,"$plot",trim(gr_plotfile)) ) 
          call execute_command_line(trim(gr_string_replace(action,"$plot",trim(gr_plotfile)) ) )
         endif
      endif

  end subroutine grende 


!== gr_color manger ========================================================================================= 
  function gr_color_manager(rgb, newindex) result(colorindex)   ! check for present index ?
                                                                ! caution multiple calls with same color
                                                                ! increments index nevertheless !
    implicit none
    double precision, intent(in)  :: rgb(3)
    integer, intent(in), optional :: newindex
    integer                       :: colorindex

    if(present(newindex)) then
      next_color_index = newindex
    else
      next_color_index = min(MAX_COLOR_INDEX,next_color_index+1) 
    endif

    if(next_color_index == MAX_COLOR_INDEX) write(*,*)"WARNING: gr_color_manager RANGE OF COLOR EXHAUSTED!"
    colorindex = next_color_index

    call gr_setcolorrep(colorindex, rgb(1), rgb(2), rgb(3))

  end function gr_color_manager

! TBD set colorscale/range

!== gr_line ================================================================================================ 
  subroutine grline( x, y, n, color, thickness, typ )
    implicit none
    double precision, intent(in)              :: x(:)
    double precision, intent(in)              :: y(:)
    integer         , intent(in), optional    :: n
    integer         , intent(in), optional    :: color
    double precision, intent(in), optional    :: thickness
    integer         , intent(in), optional    :: typ

    integer :: np

    if(present(color    )) call gr_setlinecolorind(color)
    if(present(thickness)) call gr_setlinewidth   (thickness)
    if(present(typ))       call gr_setlinetype    (typ)

    if(present(n)        ) then
      np = min(size(x), n)
    else
      np = size(x)
    endif
 
    call gr_polyline(np, x, y)

  end subroutine grline


!== grmarker ================================================================================================ 
  subroutine grsymbol( x, y, n, xerror, yerror, color, symbolsize, typ )
    implicit none
    double precision, intent(in)              :: x(:)
    double precision, intent(in)              :: y(:)
    integer         , intent(in), optional    :: n
    double precision, intent(in), optional    :: xerror(:)
    double precision, intent(in), optional    :: yerror(:)
    integer         , intent(in), optional    :: color
    double precision, intent(in), optional    :: symbolsize
    integer         , intent(in), optional    :: typ

    double precision, allocatable             :: xp(:), xm(:)

    integer :: np, nex, ney

    if(present(color     )) then
                            call gr_setmarkercolorind(color)
                            call gr_setlinecolorind(color)
    endif
    if(present(symbolsize)) then
                            call gr_setmarkersize    (symbolsize)
                            call gr_setlinewidth     (symbolsize/1.5)
    endif
    if(present(typ))        call gr_setmarkertype    (typ)

    if(present(n)        ) then
      np = min(size(x),size(y), n)    
    else
      np = min(size(x),size(y))
    endif
 
    call gr_polymarker(np, x, y)  


    if(present(xerror)) then
     nex = min(np, size(xerror))
     if(maxval(abs(xerror(1:nex))) > 0d0) then 
      allocate(xm(nex))
      allocate(xp(nex))
      xm = x(1:nex)-xerror(1:nex)
      xp = x(1:nex)+xerror(1:nex)
      call gr_herrorbars(nex, x, y, xm, xp)   
      deallocate(xm)
      deallocate(xp)
     endif
    endif

    if(present(yerror)) then
     ney = min(np, size(yerror))
     if(maxval(abs(yerror(1:ney))) > 0d0) then 
      allocate(xm(ney))
      allocate(xp(ney))
      xm = y(1:ney)-yerror(1:ney)
      xp = y(1:ney)+yerror(1:ney)
      call gr_verrorbars(ney, x, y, xm, xp)
      deallocate(xm)
      deallocate(xp)
     endif
    endif

  end subroutine grsymbol

!== grfill ================================================================================================= 
  subroutine grfill( x, y, n, color, typ )
    implicit none
    double precision, intent(in)              :: x(:)
    double precision, intent(in)              :: y(:)
    integer         , intent(in), optional    :: n
    integer         , intent(in), optional    :: color
    integer         , intent(in), optional    :: typ

    integer :: np

    if(present(color    )) call gr_setfillcolorind( color)
    if(present(color    )) call gr_setlinecolorind( color)
    if(present(typ      )) call gr_setfillintstyle (typ)

    if(present(n)        ) then
      np = min(size(x), size(y), n)
    else
      np = size(x)
    endif
 
    call gr_fillarea(np, x, y)

  end subroutine grfill


!== graxes (automatic) ===================================================================================== 
  subroutine graxes(xlabel, ylabel, title, color, option, tick_scale, text_scale)
    implicit none
    character(len=*), intent(in), optional :: xlabel
    character(len=*), intent(in), optional :: ylabel
    character(len=*), intent(in), optional :: title
    integer         , intent(in), optional :: color
    integer         , intent(in), optional :: option ! OPTION_X_LOG, OPTION_Y_LOG, OPTION_Z_LOG, OPTION_FLIP_X, OPTION_FLIP_Z
    double precision, intent(in), optional :: tick_scale
    double precision, intent(in), optional :: text_scale
!
   double precision  :: x_tick, y_tick
   double precision  :: x_org , y_org

   integer         , save  :: major_x = 1
   integer         , save  :: major_y = 1
   double precision, parameter  :: tick_size0 = 0.004d0
   double precision, save  :: tick_size= tick_size0 
   double precision, save  :: xmin, xmax, ymin, ymax
   double precision, parameter  :: textsize0 = 0.018d0 ! 0.024d0
   double precision, save  :: text_size = textsize0

   double precision        :: xlabel_x, xlabel_y
   double precision        :: ylabel_x, ylabel_y
   double precision        :: tlabel_x, tlabel_y

   integer           :: act_opt

   if(present(color)) then
     call gr_setlinecolorind(color)
     call gr_settextcolorind(color)
   endif
    
   if(present(option))      call gr_setscale(option)  
   if(present(tick_scale))  tick_size = tick_size0 * tick_scale
   if(present(text_scale))  text_size = textsize0 * text_scale

   call gr_inqwindow(xmin, xmax, ymin, ymax)

   x_org = xmin
   y_org = ymin

   call gr_inqscale(act_opt)
 
   x_tick = get_tick(xmin,xmax)  
   y_tick = get_tick(ymin,ymax) 
 
     
   call gr_setcharheight(text_size) 
   call gr_settextalign(  TEXT_HALIGN_CENTER , TEXT_VALIGN_HALF  ) 
   call gr_setlinewidth(AXIS_LINEWIDTH)
  
   call gr_axes(x_tick,y_tick,x_org,y_org, major_x,major_y, tick_size)

  ! write labels
   call gr_selntran(0)

!write(*,*)"Test axis xl:", trim(xlabel)
!write(*,*)"Test axis yl:", trim(ylabel)
  
   xlabel_x = xmin+(xmax-xmin)*0.5d0
   xlabel_y = ymin-(ymax-ymin)*XLEG_DISTANCE

   ylabel_x = xmin-(xmax-xmin)*YLEG_DISTANCE
   ylabel_y = ymin+(ymax-ymin)*0.5d0

   tlabel_x = xmin+(xmax-xmin)*TLEG_DISTANCE_X   ! to be made individual for Tit
   tlabel_y = ymax+(ymax-ymin)*TLEG_DISTANCE_Y

   if(act_opt == OPTION_X_LOG .or. act_opt == OPTION_XY_LOG) then
        xlabel_x = 10d0**(log10(xmin)+(log10(xmax)-log10(xmin))*0.5d0)
        ylabel_x = 10d0**(log10(xmin)-(log10(xmax)-log10(xmin))*YLEG_DISTANCE)
        tlabel_x = 10d0**(log10(xmin)+(log10(xmax)-log10(xmin))*TLEG_DISTANCE_X) ! to be made individual for Tit
   endif

   if(act_opt == OPTION_Y_LOG .or. act_opt == OPTION_XY_LOG) then
        xlabel_y = 10d0**(log10(ymin)-(log10(ymax)-log10(ymin))*XLEG_DISTANCE)
        ylabel_y = 10d0**(log10(ymin)+(log10(ymax)-log10(ymin))*0.5d0)
        tlabel_y = 10d0**(log10(ymax)+(log10(ymax)-log10(ymin))*TLEG_DISTANCE_Y)
   endif


   if(present(ylabel)) then
     call gr_settextpath (  TEXT_PATH_UP    )
     call gr_setcharup   (  -1d0, 0d0       )
!    call gr_textext(textsize+0.005d0 ,0.5d0,trim(grtex_filter(ylabel))//czero)
       call grtext(ylabel_x,ylabel_y,trim(ylabel)//czero)
   endif

   if(present(xlabel)) then
     call gr_settextpath (  TEXT_PATH_RIGHT )
     call gr_setcharup   (  0d0, 1d0        )
!     call gr_textext(0.5d0 ,textsize+0.005d0,trim(grtex_filter(xlabel))//czero)
      call grtext(xlabel_x,xlabel_y, trim(xlabel))
     endif

   if(present(title)) then
     call gr_settextalign(  TEXT_HALIGN_LEFT , TEXT_VALIGN_HALF  ) 
     call gr_settextpath (  TEXT_PATH_RIGHT )
     call gr_setcharup   (  0d0, 1d0        )
!     call gr_textext(0.1d0,1-2*textsize,trim(grtex_filter(title))//czero)
     call grtext(tlabel_x, tlabel_y,trim(title)//czero)
   endif


   call gr_settextalign(  TEXT_HALIGN_LEFT , TEXT_VALIGN_NORMAL  ) 

   call gr_selntran(1)
     
  end subroutine graxes

  subroutine graxes2(xlabel, ylabel, title, color, option, tick_scale, text_scale)  ! axis origin upper right corner
    implicit none
    character(len=*), intent(in), optional :: xlabel
    character(len=*), intent(in), optional :: ylabel
    character(len=*), intent(in), optional :: title
    integer         , intent(in), optional :: color
    integer         , intent(in), optional :: option ! OPTION_X_LOG, OPTION_Y_LOG, OPTION_Z_LOG, OPTION_FLIP_X, OPTION_FLIP_Z
    double precision, intent(in), optional :: tick_scale
    double precision, intent(in), optional :: text_scale

!
   double precision  :: x_tick, y_tick
   double precision  :: x_org , y_org

   integer          ,save :: major_x = 1
   integer          ,save :: major_y = 1
   double precision, parameter  :: tick_size0 = 0.004d0
   double precision, save  :: tick_size= tick_size0 
   double precision, save  :: xmin, xmax, ymin, ymax
   double precision, parameter  :: textsize0 = 0.018d0 ! 0.024d0
   double precision, save  :: text_size = textsize0

   double precision        :: xlabel_x, xlabel_y
   double precision        :: ylabel_x, ylabel_y
   double precision        :: tlabel_x, tlabel_y


   integer           :: act_opt

   if(present(color)) then
     call gr_setlinecolorind(color)
     call gr_settextcolorind(color)
   endif
   
   if(present(option)) call gr_setscale(option)  
   if(present(tick_scale))  tick_size = tick_size0 * tick_scale
   if(present(text_scale))  text_size = textsize0 * text_scale


   call gr_inqwindow(xmin, xmax, ymin, ymax)

   x_org = xmax
   y_org = ymax

   call gr_inqscale(act_opt)
 
   x_tick = get_tick(xmin,xmax)  
   y_tick = get_tick(ymin,ymax) 

   
   call gr_setcharheight(text_size) 
   call gr_settextalign(  TEXT_HALIGN_CENTER , TEXT_VALIGN_HALF  ) 

   call gr_setlinewidth(AXIS_LINEWIDTH)   
   call gr_axes(x_tick,y_tick,x_org,y_org, major_x,major_y, tick_size)

  ! write labels
   call gr_selntran(0)

   xlabel_x = xmin+(xmax-xmin)*0.5d0
   xlabel_y = ymax+(ymax-ymin)*XLEG_DISTANCE

   ylabel_x = xmax+(xmax-xmin)*YLEG_DISTANCE
   ylabel_y = ymin+(ymax-ymin)*0.5d0

   tlabel_x = xmax-(xmax-xmin)*TLEG_DISTANCE_X
   tlabel_y = ymax+(ymax-ymin)*TLEG_DISTANCE_Y

   if(act_opt == OPTION_X_LOG .or. act_opt == OPTION_XY_LOG) then
        xlabel_x = 10d0**(log10(xmin)+(log10(xmax)-log10(xmin))*0.5d0)
        ylabel_x = 10d0**(log10(xmax)+(log10(xmax)-log10(xmin))*YLEG_DISTANCE)
        tlabel_x = 10d0**(log10(xmax)-(log10(xmax)-log10(xmin))*TLEG_DISTANCE_X)
   endif

   if(act_opt == OPTION_Y_LOG .or. act_opt == OPTION_XY_LOG) then
        xlabel_y = 10d0**(log10(ymax)+(log10(ymax)-log10(ymin))*XLEG_DISTANCE)
        ylabel_y = 10d0**(log10(ymin)+(log10(ymax)-log10(ymin))*0.5d0)
        tlabel_y = 10d0**(log10(ymax)+(log10(ymax)-log10(ymin))*TLEG_DISTANCE_Y)
   endif




   if(present(ylabel)) then
     call gr_settextpath (  TEXT_PATH_DOWN    )
     call gr_setcharup   (  1d0, 0d0       )
!     call gr_textext(1-2*textsize+0.005d0 ,0.5d0,trim(grtex_filter(ylabel))//czero)
     call grtext(ylabel_x,ylabel_y,trim(ylabel)//czero)
   endif

   if(present(xlabel)) then
     call gr_settextpath (  TEXT_PATH_RIGHT )
     call gr_setcharup   (  0d0, 1d0        )
!     call gr_textext(0.5d0 ,1-3*textsize+0.005d0,trim(grtex_filter(xlabel))//czero)
      call grtext(xlabel_x, xlabel_y, xlabel)
   endif

   if(present(title)) then
     call gr_settextalign(  TEXT_HALIGN_RIGHT , TEXT_VALIGN_HALF  ) 
     call gr_settextpath (  TEXT_PATH_RIGHT )
     call gr_setcharup   (  0d0, 1d0        )
!     call gr_textext(0.9d0,1-2*textsize,trim(grtex_filter(title))//czero)
     call grtext(tlabel_x,tlabel_y,trim(title)//czero)
   endif


   call gr_settextalign(  TEXT_HALIGN_LEFT , TEXT_VALIGN_NORMAL  ) 

   call gr_selntran(1)
     
  end subroutine graxes2

  function get_tick(amin, amax) result(tick)
    implicit none
    double precision, intent(in) :: amin 
    double precision, intent(in) :: amax
    double precision  :: tick
  
    integer  :: m, u
  
    m      =   floor(log10(amax-amin))
    tick   =   10.0d0**m
    u      =   (amax-amin)/tick


    select case(u)
      case(5:10) 
        tick = tick*2
      case(3:4) 
        tick = tick
      case(2)
        tick = tick/2
      case(0:1)
        tick = tick/4
      case default
        tick = tick
    end select
   
  end function get_tick

!== gr_text ================================================================================================ 
  
  subroutine grtext(xwc, ywc, txt, colorindex, xtab)
   implicit none
   double precision,   intent(in)           :: xwc, ywc  ! world coordinates for text
   character(len=*),   intent(in)           :: txt
   integer         ,   intent(in), optional :: colorindex
   double precision,   intent(in), optional :: xtab

   integer, external :: gr_wctondc ! world coordinate to normalized device coordinate
   integer           :: i,j   
   double precision  :: xt , yt

   if(present(colorindex))  call gr_settextcolorind( colorindex )

   xt = xwc
   yt = ywc


   i = gr_wctondc(xt,yt)
! experimental use "=" as tabulator (only one) to be improved TBD...
   if(present(xtab) .and. index(txt,"=") > 0) then
       j = index(txt,"=")
       call gr_text(xt,yt,txt(1:j-1)//czero)
       xt = xwc + xtab
       yt = ywc
       i = gr_wctondc(xt,yt)
       call gr_text(xt,yt,trim(txt(j:))//czero)
   else
       if(txt(1:1) .ne. "$") then
         call gr_text(xt,yt,trim(txt)//czero)
       else
          if(txt(2:2) .ne. "$") then
            call gr_textext(xt,yt,trim(grtex_filter(txt(2:)))//czero)   ! one $ at start
          else
            call gr_mathtex(xt,yt,trim(grtex_filter(txt(3:)))//czero)   ! two $$ at start
          endif
       endif
   endif
    
  end subroutine grtext

!-------------------------------------------------------------------------------------------------- 
     
  function grtex_filter(sinitial) result(sfinal)
    implicit none
    character(len=*), intent(in)   :: sinitial
    character(len=2*len(sinitial)) :: sfinal

    sfinal = gr_string_replace(sinitial,"_","\_")
    sfinal = gr_string_replace(sfinal  ,"$","\$x")
!    sfinal = gr_string_replace(sfinal  ,"%","\%")
    sfinal = gr_string_replace(sfinal  ,"#","\#")
!    sfinal = gr_string_replace(sfinal  ,">","$>$")
!    sfinal = gr_string_replace(sfinal  ,"<","$<$")
    sfinal = gr_string_replace(sfinal  ,"/","\/")

  end function grtex_filter


     
  function grtitle_filter(sinitial) result(sfinal)
    implicit none
    character(len=*), intent(in)   :: sinitial
    character(len=2*len(sinitial)) :: sfinal

    integer :: i
    sfinal = " "
    do i=1,len_trim(sinitial)
     sfinal(i:i) = sinitial(i:i)
     if(sinitial(i:i) == " ") sfinal(i:i) = "_"
     if(sinitial(i:i) == "(") sfinal(i:i) = "["
     if(sinitial(i:i) == ")") sfinal(i:i) = "]"
     if(sinitial(i:i) == "$") sfinal(i:i) = "t"
     if(sinitial(i:i) == "\") sfinal(i:i) = "_"
    enddo
   

  end function grtitle_filter


!== gr_fmt ==================================================================================================  
   
  !function gr_fmt_int(num,digit, digit2) result(sfinal) - paz - unused digit2
  function gr_fmt_int(num,digit) result(sfinal) 
    implicit none
    integer, intent(in)            :: num
    integer, intent(in)            :: digit
    !integer, intent(in), optional  :: digit2 - unused parameter


    character(len=16)              :: sfinal
    character(len=16)              :: ifmt

    write(ifmt,'("(i",i0,")")') digit
    write(sfinal,ifmt)     num

  end function gr_fmt_int

!-------------------------------------------------------------------------------------------------- 
   
  function gr_fmt_dbl(num,digit, digit2) result(sfinal)
    implicit none
    double precision, intent(in)   :: num
    integer, intent(in)            :: digit
    integer, intent(in), optional  :: digit2

    character(len=16)              :: sfinal
    character(len=16)              :: ifmt

    if(present(digit2)) then
       write(ifmt,'("(f",i0,".",i0,")")') digit, digit2
    else
       write(ifmt,'("(f",i0,".",i0,")")') max(3,digit), max(3,digit)/2
    endif
    write(sfinal,ifmt)     num

  end function gr_fmt_dbl

!== gr_fmt ==================================================================================================  

  function gr_string_replace(string, substring1, substring2) result(new_string)
    character(len=*), intent(in)         :: string
    character(len=*), intent(in)         :: substring1
    character(len=*), intent(in)         :: substring2
    character(len=len(string)+(len(string)/len(substring1)+1)*len(substring2)) :: new_string
    integer  :: i, k, l, i1, i2, len0, len1, len2

    l    = len(string)
    len0 = len_trim(string)
    len1 = len_trim(substring1)
    len2 = len_trim(substring2)

    new_string = ' '
    i1 = 1
    i2 = l

    do k=1,l
        i = index(string(i1:l),substring1(1:len1))
        if(i == 0) then
            new_string = trim(new_string)//string(i1:l)
            exit
        endif

        i = i + i1 - 1
        new_string = trim(new_string)//string(i1:i-1)//substring2(1:len2)
        i1 = i+len1
    enddo

    new_string = trim(new_string)
    return
  end function gr_string_replace



 END  module gr_interface 
!example!  
!example! 
!example! 
!example! 
!example!                                    
!example! !                                                                       
!example! !   gfortran grtestmod.f90 -Wl,-rpath,/usr/local/gr/lib  -L/usr/local/gr/lib -lGR
!example! !                                                                       
!example! program grtestmod 
!example! !
!example!  use gr_interface
!example! 
!example! implicit none
!example! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!example!       integer          :: ici
!example!       integer          :: red, blue, green, yellow, black
!example! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!example! !                                                                       
!example!       double precision xd(100), yd(100), zd(100) 
!example!       double precision x(200), y(200), z(200, 200) 
!example!       double precision h(20) 
!example!       integer i , ind(40)
!example! 
!example!       call grstart(plotfile="test2.pdf",boundary=0.15d0)
!example! !     call grstart(boundary=0.15d0)
!example! 
!example!       call gr_setwindow(0d0, 20.0D0, -5.0d0, 5.0D0) 
!example! 
!example! ! example for spefically defined colors from the rgb-listing 
!example! ! (simple colors are predefined: GR_RED, GR_BLUE, ....)
!example!       black  = gr_color_manager(RGB_BLACK)
!example!       red    = gr_color_manager(RGB_RED)
!example!       green  = gr_color_manager(RGB_GREEN)
!example!       blue   = gr_color_manager(RGB_BLUE)
!example!       yellow = gr_color_manager(RGB_YELLOW)
!example! !  ... etc for all desired colors
!example! 
!example!       call graxes ("xlabel ist hier","und ylabel",'Ueberschrift',GR_BLUE,OPTION_LINEAR)
!example!       call graxes2("xlabel2 ist hier","und ylabel2",'Ueberschrift2',red,OPTION_LINEAR)
!example! 
!example! 
!example! 
!example!       x = [ (0.1d0*i, i=1,200) ]
!example!       y = 4*sin(3*(x-4)) 
!example! 
!example!       call gr_polyline(100,x(1:100),y(1:100)) ! this is the generic ployline
!example! 
!example!       call grline(x=x, y=y, n=100, color=yellow, thickness=7d0,   typ=LINETYPE_DASHED)
!example!       call grline(x=x, y=y, n=2,   color=blue,   thickness=0.5d0, typ=LINETYPE_SOLID)
!example! 
!example!       call grline(x,y)     ! simplest call, now using all 200 elements of x,y
!example! 
!example! ! == draw a symbol every 5th point of x,y
!example!       ind = [(5*i,i=1,ubound(x,dim=1)/5)]
!example!       call grsymbol(x=x(ind), y=y(ind), xerror    = y(ind)*0.2d0, yerror = y(ind)*0.1d0, &
!example!                                         symbolsize=1.5d0,           color=GR_MAGENTA)
!example! ! == and a  filled rectangle
!example! !      call grfill(x=2*[0d0,1d0,1d0,0d0,0d0]+3,y=[0d0,0d0,1d0,1d0,0d0]*2,color=green,typ=SOLID)
!example!                                                                       
!example! !     call gr_updatews() 
!example! !     call grende   ! simplest for of use
!example!       call grende("open $plot ; cp $plot lastplot.pdf")  ! $plot is symbolic for then name of the 
!example!                                                          ! created plotfile
!example!                         
!example! end program grtestmod          

       subroutine splot (doplo)
!      ================  scan-plotting
       use gr_interface

       use new_com
       ! use cincoc
       use cdata
!       use outlev
       use theory
       use selist
       use fslist
       use therrc
       use thparc
       use constants
       use formul
       implicit none

       character(len=40)   :: LAST_PPAR = "last_plotsetting"

       double precision    :: e(size(xwerte(:,1)))
       double precision    :: x(size(xwerte(:,1)))
       double precision    :: y(size(xwerte(:,1)))
       integer,save :: icolo(size(inpar)) = 0
       integer,save :: irecv(size(inpar))
       integer :: i
       integer,save :: isymb(1:size(inpar)) = [(i,i=1,size(inpar))]+1   !! check Markertypes 
       integer,save :: ifrec(1:size(inpar))

       double precision, save :: linewidth_scaling(1:size(inpar)) = 1d0
       double precision, save :: sysize_scaling   (1:size(inpar)) = 1d0
       integer         , save :: linetype(1:size(inpar)) = LINETYPE_SOLID


       double precision, save       :: p0_scale(size(inpar)) = 1.0
       double precision             :: p_scale(size(inpar))

       integer          :: inew

       logical :: ptex=.true.
       logical :: paplo=.true.
       logical :: doplo
       logical :: fitplo=.true.
       logical :: errplo=.true.
       logical :: paxis=.true.
       logical :: taxis=.true.

       logical :: log_x=.false.
       logical :: log_y=.false.
!      --- doplo = false  means: set parameters only ---
       character(len=80) :: xtext,ytext
       character(len=80) :: buf
       character(len=12) :: tag, stunde
       character(len=12) :: tx,sx

      double precision, save :: xmin=0.d0, xmax=1.d0, ymin=0.d0 ,ymax=1.d0
      integer         , save :: nkurv=0
      double precision, save :: txsize=0.35d0
      double precision, save :: sysize=0.6d0
      double precision, save :: fyskip=1.2d0
      double precision, save :: epscl  = 0.001
      integer         , save :: icol0  = 1
      double precision, save :: txsizt = 0.23d0 
      double precision, save :: xtshft = 0.0d0
      double precision, save :: ytshft = 0.0d0
    
      double precision :: thcline_thickness = 0.5d0
      double precision :: datline_thickness = 0.75d0

      double precision, save :: ax_text_scale = 1d0
      double precision, save :: ax_tick_scale = 1d0

      double precision :: yhigh, ylow, xskip

      double precision ::  ytxs, yma_s, ymi_s, xtxs, ytx
      double precision ::  yepl, yeml, xtx, xmi_s,xma_s
      integer irfcu, j, icco, ik, ip, ircu, nsy, npic, nnpi, it
      integer npicf, npar, nfkurv, nco, lxx, lyy, l, ith, ircf

      integer, save :: iplevel = 0
      integer       :: axis_option

      integer       ::  ifont = 0
      integer       ::  ubild
      logical       ::  fileda
      
      logical, save :: initial = .true.



      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= plot                                                                       ='
       write(6,*)'=    plots selected records                                                  ='
       write(6,*)'=    parameters:                                                             ='
       write(6,*)'=      xmin  <val>    :  start of x-plotting range                           ='
       write(6,*)'=      xmax  <val>    :  end   of x-plotting range                           ='
       write(6,*)'=      ymin  <val>    :  start of y-plotting range                           ='
       write(6,*)'=      ymax  <val>    :  end   of y-plotting range                           ='
       write(6,*)'=      symb  <list>   :  integer synbmbol list (sequence = selected)         ='
       write(6,*)'=      sysize <val>   :  sets symbol size                                    ='
       write(6,*)'=      icolo <list>   :  integer color    list (sequence = selected)         ='
       write(6,*)'=      ltype <list>   :  linetypes                                           ='
       write(6,*)'=      lwid  <list>   :  linewidths                                          ='
       write(6,*)'=      sywid <list>   :  symbol size scalings                                ='
       write(6,*)'=      errors         :  adds error bars                                     ='
       write(6,*)'=      noerrors       :  suppress error bars                                 ='
       write(6,*)'=      parplo         :  sets parameter value listing                        ='
       write(6,*)'=      parlev         :  sets level for (more) parameter lsiting             ='
       write(6,*)'=      noparplo p1 p2.:  suppress parameter listing except for p1 p2 ...     ='
       write(6,*)'=      txsize <val>   :  set textsize (legend)                               ='
       write(6,*)'=      axtxsize <val> :  scale axis script text size                         ='
       write(6,*)'=      axticlen <val> :  scale axis tick size (neg=outbound)                 ='
       write(6,*)'=      xlegdist <val> :  distance of x-axis name from axis                   ='
       write(6,*)'=      ylegdist <val> :  distance of y-axis name from axis                   ='
       write(6,*)'=      tit_x <val>    :  distance of title from axis left                    ='
       write(6,*)'=      tit_x <val>    :  distance of title from axis top                     ='
       write(6,*)'=      flinewd  <val> :  fit linewidth                                       ='
       write(6,*)'=      dfinewd  <val> :  data linewidth                                      ='
       write(6,*)'=      lin_x | log_x  :  lin or log x scale                                  ='
       write(6,*)'=      lin_y | log_y  :  lin or log y scale                                  ='
       write(6,*)'=      def            :  at first run start with default                     ='
!       write(6,*)'=      last           :  read parameters form last session                   ='
       write(6,*)'=      # <al>         :  auto num picture store  ON start at <val>, neg=OFF  ='
       write(6,*)'=        for further options see manual ....                                 ='
       write(6,*)'= HINTS:  (prior to plot)                                                    ='
       write(6,*)'=       use the:  title    command to set a plot title                       ='
       write(6,*)'=       use the:  rename   command to change x-axis and y-axis names         ='
       write(6,*)'=                 rename   names starting with $ may contain some tex codeing='
       write(6,*)'= PRINTING/SAVING:                                                           ='
       write(6,*)'=       if there is some title (use tit abc...) then                         ='
       write(6,*)'=       plot are stored to dtrplot_<title>-#.pdf                             ='
       write(6,*)'=       (blanks etc. are replaced by _ etc.) # is the plotnumber shown also  ='
       write(6,*)'=       on the lower left cormer of the plot                                 ='
       write(6,*)'=       and always copied to last_datreat_plot.pdf                           ='
       write(6,*)'=       to stop creation of auto numbered series datreat_plot##.pdf  use plot # <neg>'
       write(6,*)'=       to restart creation of auto numbered series datreat_plot##.pdf  use plot # <n>'
       write(6,*)'= VERSION 3.0                                                                ='
       write(6,*)'=============================================================================='
       return
      endif


! ----- parameter retrieving from stack -----
      nkurv  = 0
      nfkurv = 0
      if(inames.eq.0) then
!                     ----> assume that a list of number numors
       if(ipars.gt.0) then
          nkurv = 0
          nsel  = 0
          if(ipars.gt.size(inpar)) ipars = size(inpar)
          do 2 i=1,ipars
           irecv(i) = rpar(i) + 0.0001
    2     continue
          nkurv = ipars
        endif
       else
!      -----> decode by names
!write(*,*)"t00:", inames
        do 3 i=1,inames
          j = inapa(i)
          if(vname(i).eq.'xmin    ') xmin = rpar(j)
          if(vname(i).eq.'xmax    ') xmax = rpar(j)
          if(vname(i).eq.'ymin    ') ymin = rpar(j)
          if(vname(i).eq.'ymax    ') ymax = rpar(j)
          if(vname(i).eq.'text    ') ptex = .true.
          if(vname(i).eq.'txon    ') ptex = .true.
          if(vname(i).eq.'notext  ') ptex = .false.
!   ---> txon, txoff are from old version and not documented any more <-
          if(vname(i).eq.'txoff   ') ptex  = .false.
          if(vname(i).eq.'epscl   ') epscl = rpar(j)
          if(vname(i).eq.'txsize  ') txsize= rpar(j)
          if(vname(i).eq.'legsize ') txsizt= rpar(j)
          if(vname(i).eq.'legx    ') xtshft= rpar(j)
          if(vname(i).eq.'legy    ') ytshft= rpar(j)
          if(vname(i).eq.'color   ') icol0 = rpar(j) +0.0001
          if(vname(i).eq.'sysize  ') sysize= rpar(j)
          if(vname(i).eq.'parlev  ') iplevel = Nint(rpar(j))
          if(vname(i).eq.'parplo  ') paplo=.true.
          if(vname(i).eq.'noparplo') paplo=.false.
          if(vname(i).eq.'errplo  ') errplo=.true.
          if(vname(i).eq.'noerrplo') errplo=.false.
          if(vname(i).eq.'errors  ') errplo=.true.
          if(vname(i).eq.'noerrors') errplo=.false.
          if(vname(i).eq.'axis    ') paxis =.true.
          if(vname(i).eq.'noaxis  ') paxis =.false.
          if(vname(i).eq.'txaxis  ') taxis =.true.
          if(vname(i).eq.'notxaxis') taxis =.false.
          if(vname(i).eq.'fits    ') fitplo = .true.
          if(vname(i).eq.'nofits  ') fitplo = .false.
          if(vname(i).eq.'log_x   ') log_x  = .true.
          if(vname(i).eq.'log_y   ') log_y  = .true.
          if(vname(i).eq.'lin_x   ') log_x  = .false.
          if(vname(i).eq.'lin_y   ') log_y  = .false.
          if(vname(i).eq.'def     ') initial= .false.
          if(vname(i).eq.'last    ') initial= .true.
!
          ibild               = get_named_value("#          ",ibild,inew)
          ifont               = get_named_value("font       ",ifont,inew)
          ax_text_scale       = get_named_value("axtxsize   ",ax_text_scale,inew)
          ax_tick_scale       = get_named_value("axticlen   ",ax_tick_scale,inew)

          thcline_thickness   = get_named_value("flinewd   ",thcline_thickness,inew)
          datline_thickness   = get_named_value("dlinewd   ",datline_thickness,inew)


! experimental:
          XLEG_DISTANCE   = get_named_value("xlegdist   ",XLEG_DISTANCE,inew)
          YLEG_DISTANCE   = get_named_value("ylegdist   ",YLEG_DISTANCE,inew)
          TLEG_DISTANCE_X = get_named_value("tit_x      ",TLEG_DISTANCE_X,inew)
          TLEG_DISTANCE_Y = get_named_value("tit_y      ",TLEG_DISTANCE_Y,inew)

          
! to be modernized
          if(vname(i).eq.'ltype    ' ) then
            nsy = 0
            do  l=1,inpar(i)
             nsy   = nsy   + 1
             if(nsy.gt.size(inpar)) exit
             linetype(nsy) = nint(abs( rpar(j) ))
             j = j + 1
            enddo
          endif
           
! to be modernized
          if(vname(i).eq.'lwid    ' ) then
            nsy = 0
            do  l=1,inpar(i)
             nsy   = nsy   + 1
             if(nsy.gt.size(inpar)) exit
             linewidth_scaling(nsy) =  rpar(j)
             j = j + 1
            enddo
          endif
 

!! to be modernized
          if(vname(i).eq.'symb    ' .or. vname(i).eq.'isymb   ' ) then
            nsy = 0
            do 49 l=1,inpar(i)
             nsy   = nsy   + 1
             if(nsy.gt.size(inpar)) goto 29
             isymb(nsy) = nint(abs( rpar(j) ))
             j = j + 1
   49       continue
          endif
   29     continue
!
!! to be modernized
 fsywd:  if(vname(i).eq.'sywid    ') then
            nsy = 0
            do l=1,inpar(i)
             nsy   = nsy   + 1
             if(nsy.gt.size(inpar)) exit
             sysize_scaling(nsy) =  rpar(j)
             j = j + 1
            enddo
          endif fsywd
   
!

          if(vname(i).eq.'colo    ' .or. vname(i).eq.'icolo   ') then
            nco = 0
            do 59 l=1,inpar(i)
             nco   = nco   + 1
             if(nco.gt.size(inpar)) goto 39
             icolo(nco) = nint( abs(rpar(j)))
             j = j + 1
   59       continue
          endif
   39     continue
!
          if(vname(i).eq.'sc      ') then
            nsel = 0
            nkurv= 0
            do 5 l=1,inpar(i)
             nkurv = nkurv + 1
             if(nkurv.gt.size(inpar)) goto 31
             irecv(nkurv) = Nint(rpar(j))
             j = j + 1
    5       continue
          endif
   31     continue
!
          if(vname(i).eq.'fsc     ') then
            nfsel = 0
            nfkurv= 0
            do 641 l=1,inpar(i)
             nfkurv = nfkurv + 1
             if(nfkurv.gt.size(inpar)) goto 64
             ifrec(nfkurv) = Nint(rpar(j))
             j = j + 1
  641       continue
          endif
   64     continue
          if(nfsel.eq.0) then
            if(nfkurv.gt.0) then
               call fsrch(ifrec,nfkurv)
            endif
          endif
!
    3   continue
      endif
!
      if (.not.doplo) then
        return
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(initial) then 
        call readplotpar()
        initial = .false.
      endif
      call writeplotpar()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
      if(nsel.eq.0) then
        if(nkurv.gt.0) then
          call search(irecv,nkurv)
!       -----------------------------> select the items to be plotted
        else
          write(6,*)'no scan selected, no plot possible !'
          return
        endif
      endif
!
      if (.not.doplo) then
        return
      endif
! ---->                                             <-----------
      nkurv = nsel
! ----> ich denke, das muss so sein !?!?!?!  (m.s.) <-----------
      if (nkurv.eq.0) then
        write(6,*)'there is no curve selected => no plot !'
        return
      endif
!
scl:   if(found('scaled  ')) then
         p0_scale(1) = getval('scaled  ',dble(p0_scale(1)),inew)
         do i=2,size(inpar)
          p0_scale(i) = valnxt(dble(p0_scale(i)),inew)
          if(inew.eq.0) exit
         enddo
          p_scale(1:size(inpar)) = p0_scale(1:size(inpar))
       else
          p_scale(1:size(inpar)) = 1.0d0
       endif scl

!
!       if(ibild>0) ibild = ibild + 1
       if(ibild>=0) then         
          inquire(file="datreat_plot_counter",exist=fileda)
          if(fileda) then
            open(newunit=ubild,file="datreat_plot_counter")
            read(ubild,*) ibild
            ibild = ibild + 1
            rewind(ubild)
            write(ubild,*) ibild
            close(ubild)
          else
            open(newunit=ubild,file="datreat_plot_counter")
            ibild = ibild+1
            write(ubild,*) ibild
            close(ubild)
          endif
       endif
!
! ---- set frame & scales ----
! !-> preliminary fix       
       if(len_trim(title) == 0) title=" "

!       if(len_trim(title) > 0) then
!          call grstart("dtrplot.pdf")           !>neu, TBD use frlux ....
          write(*,'(a,i0,a)')"Plot title(#",ibild-1,"): "//"dtrplot_"//trim(grtitle_filter(title))//".pdf"
          write(buf,'(i0)')ibild-1
          call grstart("dtrplot_"//trim(grtitle_filter(title))//"-"//trim(buf)//".pdf") !>neu, TBD use frlux ....
!       else
!          call grstart 
!       endif

       call gr_setwindow(xmin,xmax,ymin,ymax)      !>neu
       call gr_settextfontprec(3, ifont )   ! ???? 
       call gr_setclip(0)

         xtext = xname(isels(1))
         xmi_s = xmin
         xma_s = xmax
         ytext = yname(isels(1))
         ymi_s = ymin
         yma_s = ymax

                                     axis_option = OPTION_LINEAR
         if(log_x .and. .not. log_y) axis_option = OPTION_X_LOG
         if(log_y .and. .not. log_x) axis_option = OPTION_Y_LOG
         if(log_x .and.       log_y) axis_option = OPTION_XY_LOG


        call gr_setwindow(xmi_s,xma_s,ymi_s,yma_s)    !> neu
!      write(*,*)' make axes .....'
       if(.not.taxis) then
         lxx = 0
         lyy = 0
       else
         lxx = len_trim(xtext)
         lyy = len_trim(ytext)
       endif

!write(*,*)"test1:",xtext,ytext,lxx,lyy
!?        if(paxis) call graxs(lopt,option,lxx,xtext,lyy,ytext)
       if(paxis) then
            call graxes (trim(xtext),trim(ytext),trim(title),GR_BLACK, axis_option, &
                         ax_tick_scale, ax_text_scale * graspekt) !> neu
            call graxes2 (" "," "," ",GR_BLACK, axis_option, &
                         -ax_tick_scale, ax_text_scale*0.001d0) !> neu
       endif

!
!
! ---- plot the selected fit-kurves ----
!
       nfkurv = nfsel
       if(fitplo) then
       do 70 i=1,nfkurv
        irfcu = isfits(i)
        npicf = nwert(irfcu)
        icco = mod(icolo(i),7) + 1
        nnpi = 0
        do 70010 j=1,npicf
          if(xwerte(j,irfcu).ge.xmin.and.xwerte(j,irfcu).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,irfcu) * p_scale(i)
!!          if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).lt.ymin) y(nnpi) = ymin
            if(y(nnpi).gt.ymax) y(nnpi) = ymax
            x(nnpi) = xwerte(j,irfcu)

          endif
70010     continue
           call grline(x=x, y=y, n=nnpi,   color=icco,   thickness=thcline_thickness, typ=LINETYPE_SOLID) 
!?         call grln(x,y,nnpi)   
   70   continue
       endif
!
!
! ---- plot datarecords ----
!
       nkurv = nsel
       do 20 i=1,nkurv
        ircu = isels(i)
        ircf = ifits(i)
        icco = mod(icolo(i),7) + 1
        npic = nwert(ircu)
        if(fitplo) then
        if(ircf.ne.0) then
!                     ----> plot the fitted data automatically
!                           this is after a fit-command until you
!                           select new curves
        npicf = nwert(ircf)
        nnpi  = 0
        do 20010 j=1,npicf
          if(xwerte(j,ircf).ge.xmin.and.xwerte(j,ircf).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,ircf) * p_scale(i)
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).gt.ymax) y(nnpi) = ymax+(ymax-ymin)*0.02
            x(nnpi) = xwerte(j,ircf)
          endif
20010     continue
          call grline(x=x, y=y, n=nnpi,   color=icco,   thickness=thcline_thickness, typ=LINETYPE_SOLID)

        endif
       endif
!
!                  ---> prepare data to plot
        nnpi = 0
        do 2001 j=1,npic
          if(xwerte(j,ircu).ge.xmin.and.xwerte(j,ircu).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,ircu) * p_scale(i)
            x(nnpi) = xwerte(j,ircu) 
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).gt.ymax) y(nnpi) = ymax+(ymax-ymin)*0.02
            e(nnpi) = yerror(j,ircu) * p_scale(i)


          endif
 2001   continue
!
!
!
! --- plot ---
!       if (numor(ircu).gt.0) then
!           icco=mod(icolo(i),7) + 1

!          ----------- plot a dataline -------
           if(isymb(i).eq.0) then
             call grline(x=x, y=y, n=nnpi,   color=icco,   thickness=datline_thickness*linewidth_scaling(i),&
                                                           typ=linetype(i)) !> neu
           else
             call grsymbol(x=x, y=y,  n=nnpi, &    ! xerror    = y(ind)*0.2d0, yerror = y(ind)*0.1d0, &
                                      color=icco, symbolsize= sysize * sysize_scaling(i), &
                                      typ= MARKERTYPE(isymb(i)) )
           endif

           if(errplo) then
! old style
             do ik=1,nnpi       
               yepl = y(ik)+e(ik)
               yeml = y(ik)-e(ik)
              call grline(x=[x(ik),x(ik)],y=[yeml,yepl],n=2,thickness=datline_thickness, typ=LINETYPE_SOLID)   
             enddo
! new
             call grsymbol(x=x, y=y,  n=nnpi, yerror = e, &
                            color=icco, symbolsize= sysize* sysize_scaling(i) , typ= MARKERTYPE(max(1,isymb(i))) )
           endif

!>             call grnwpn(1)
!       endif
!
   20  continue
!
!

! ---- textpart ----
       if(ptex) then

        call gr_settextalign(TEXT_HALIGN_LEFT , TEXT_VALIGN_BOTTOM ) 
!
        txsizt = txsize * 23d0/35d0    ! legacy factors may be eliminated in refactoring
! --- title ---
         xtx = xmin + 0.1 * (xmax-xmin)
         ytx = ymax
         if(found('scaled  ')) then
            call grtext(xtx,ytx-0.05d0*(ymax-ymin),"SCALED",GR_RED)  !> neu
         endif
! ---- set textwindow ----
         call gr_setwindow( DEVICE_LOWER_LEFT_X*100 , DEVICE_UPPER_RIGHT_X*100 , &
                            DEVICE_LOWER_LEFT_Y*100 , DEVICE_UPPER_RIGHT_Y*100 )
    
         call gr_setcharheight(txsizt/35*graspekt)   !> neu
         xtx =  DEVICE_UPPER_RIGHT_X*100 + 2.0  + xtshft
         ytx =  DEVICE_UPPER_RIGHT_Y*100 + 2.0  + ytshft

         yhigh = ytx
         ylow  = DEVICE_LOWER_LEFT_Y*100 - (DEVICE_UPPER_RIGHT_Y- DEVICE_LOWER_LEFT_Y)*DEFAULT_WC_MARGIN*100
         xskip = 42 * fyskip * txsizt 
!
! - plot date and time info:
!
       call date_and_time(tag,stunde)

       tx = tag(7:8)//'-'//tag(5:6)//'-'//tag(1:4)
       sx = stunde(1:2)//':'//stunde(3:4)//':'//stunde(5:6)
!       xtext = tx//'  '//sx
       write(xtext,'(a,i0)') tx//'  '//sx//' #',ibild-1
       call grtext(-DEFAULT_WC_MARGIN*150*fyskip * txsizt ,ylow,trim(xtext),GR_BLACK) !> neu

! ---- plot theory parameters ----
         if(ntheos.ne.0) then
           do 115 it = 1,ntheos
             ith = nthtab(it)
             if(multflg(it).eq.1) then
                write(xtext,'(8htheory* ,a8)')thenam(ith)
             else
                write(xtext,'(8htheory+ ,a8)')thenam(ith)
             endif
             if(thenam(ith)=="eval    ") write(xtext,'(a,":",a)')trim(xtext),trim(yfitform)
             call grtext(xtx,ytx,trim(xtext),GR_BLACK) !> neu

             call advance_text(xtx,ytx,yhigh,ylow,xskip)

             if(thrapar(it).ne.'        ') then
               write(xtext,'(4hfor ,a8,2f12.6)')thrapar(it),thramin(it),thramax(it)
               call grtext(xtx,ytx,trim(xtext),GR_BLACK,8*txsize) !> neu
              call advance_text(xtx,ytx,yhigh,ylow,xskip)

             endif
             npar = nthpar(ith)
             if(npar.ne.0) then
               do 117 ip = 1,npar
               write(xtext,'(a8,1h=,1es12.5,2h+-,es9.2,es8.1)')thparn(ip,ith),thparx(ip,it),therro(ip,it),thpsca(ip,it)
               call grtext(xtx,ytx,trim(xtext),GR_BLACK,8*txsize) !> neu
               call advance_text(xtx,ytx,yhigh,ylow,xskip)

  117          continue
            endif
  115     continue
         endif
! ---- plotted items ----
         do 101 i=1,nkurv
           ircu = isels(i)
           write(xtext,'(a8,i14,a7,es13.6)') name(ircu),numor(ircu),' scale ',p_scale(i)
           xtxs = xtx - 2*txsizt
           ytxs = ytx + txsizt / 2
           icco=mod(icolo(i),7) + 1
           if(isymb(i).ne.0.and.numor(ircu).gt.0) then
               call grsymbol(x=[xtxs,xtxs], y=[ytxs,ytxs], n=2, &    ! xerror    = y(ind)*0.2d0, yerror = y(ind)*0.1d0, &
                                     symbolsize=sysize, color=icco, typ= MARKERTYPE(isymb(i)))
           else
             call grtext(xtxs,ytxs,"-",icco) !> neu
           endif
           call grtext(xtx,ytx,trim(xtext),icco) !> neu

           call advance_text(xtx,ytx,yhigh,ylow,xskip)

           call grtext(xtx,ytx,trim(coment(ircu)),icco) !> neu

           call advance_text(xtx,ytx,yhigh,ylow,xskip)
           if(paplo) then
           do 1012 l=1,nopar(ircu)
            if(params_display_level(l,ircu) > iplevel) cycle
           xtext = " "
           write(xtext,'(a8,1h=,1es14.6)')napar(l,ircu),params(l,ircu)
           call grtext(xtx,ytx,trim(xtext),icco,10*txsizt) !> neu
           call advance_text(xtx,ytx,yhigh,ylow,xskip)
 1012      continue
           else
 dap01:     do l=1,nopar(ircu)
             if(params_display_level(l,ircu) > iplevel) cycle dap01
             if(found(napar(l,ircu)//' ')) then
              xtext = " "
              write(xtext,'(a8,1h=,1es14.6)')napar(l,ircu),                   &
     &                                  params(l,ircu)
              call grtext(xtx,ytx,trim(xtext),icco,10*txsizt) !> neu
              call advance_text(xtx,ytx,yhigh,ylow,xskip)
            endif
           enddo dap01
           endif

           ircu = ifits(i)
           if(ircu.gt.0) then
             xtext = " "
             write(xtext,'(a8,i14,a7,1es13.6)') name(ircu),numor(ircu),' scale ',p_scale(i)
             xtxs = xtx - 2*txsizt
             ytxs = ytx + txsizt / 2
             icco=mod(icolo(i),7) + 1

             if(isymb(i).ne.0.and.numor(ircu).gt.0) then
               call grsymbol(x=[xtxs], y=[ytxs], n=1, &    ! xerror    = y(ind)*0.2d0, yerror = y(ind)*0.1d0, &
                                     symbolsize=sysize, color=icco,typ= MARKERTYPE(isymb(i)) )
             else
               call grtext(xtxs,ytxs,"-",icco) !> neu
             endif

             call grtext(xtx,ytx,trim(xtext),icco) !> neu
             call advance_text(xtx,ytx,yhigh,ylow,xskip)
             call grtext(xtx,ytx,trim(coment(ircu)),icco,10*txsizt) !> neu
             ytx = ytx - fyskip  * txsizt
             if(paplo) then
 dap1:      do  l=1,nopar(ircu)
               if(params_display_level(l,ircu) > iplevel) cycle dap1
               write(xtext,'(a8,1h=,1es14.6)')napar(l,ircu),params(l,ircu)
               call grtext(xtx,ytx,trim(xtext),icco,10*txsizt) !> neu
               call advance_text(xtx,ytx,yhigh,ylow,xskip)
             enddo dap1
             else
 dap2:      do l=1,nopar(ircu)
              if(params_display_level(l,ircu) > iplevel) cycle dap2
              if(found(napar(l,ircu)//' ')) then
                write(xtext,'(a8,1h=,1es14.6)')napar(l,ircu),                 &
     &                                   params(l,ircu)

                call grtext(xtx,ytx,trim(xtext),icco,10*txsizt) !> neu
                call advance_text(xtx,ytx,yhigh,ylow,xskip)
              endif
             enddo dap2
             endif
           endif
  101    continue
      endif


!!! finalize and save pictures on disc 
      if(ibild > 0 .and. len_trim(title)==0) then
         write(buf,'("datreat_plot",i0,".pdf")')ibild-1
      else
         buf = "last_datreat_plot.pdf"
      endif
      LAST_DTR_PLOT = buf
      call grende(trim(PDF_VIEWER_COMMAND)//' $plot; cp -f $plot '//trim(buf))  
                                            ! $plot is symbolic for then name of the 
                                            ! created plotfile 
!
      CONTAINS
       subroutine advance_text(xt,yt,yhigh,ylow,xskip)
        implicit none
        double precision, intent(inout) :: xt, yt
        double precision, intent(in)    :: yhigh, ylow
        double precision, intent(in)    :: xskip
 
        if(yt - fyskip * txsizt  < ylow) then
           yt = yhigh
           xt = xt + xskip
        else
           yt = yt - fyskip * txsizt 
        endif

       end subroutine advance_text



       subroutine readplotpar()
         implicit none
         integer :: ioppar, nl
         logical :: lastppar_da
         inquire(file=LAST_PPAR,EXIST=lastppar_da)
         if(.not. lastppar_da) return
         
         open(newunit=ioppar, file=LAST_PPAR)
          read(ioppar,*,end=99,err=99)  xmin                   
          read(ioppar,*,end=99,err=99)  xmax                
          read(ioppar,*,end=99,err=99)  ymin                
          read(ioppar,*,end=99,err=99)  ymax                
          read(ioppar,*,end=99,err=99)  ptex                
          read(ioppar,*,end=99,err=99)  epscl               
          read(ioppar,*,end=99,err=99)  txsize              
          read(ioppar,*,end=99,err=99)  txsizt              
          read(ioppar,*,end=99,err=99)  xtshft              
          read(ioppar,*,end=99,err=99)  ytshft              
          read(ioppar,*,end=99,err=99)  icol0               
          read(ioppar,*,end=99,err=99)  sysize              
          read(ioppar,*,end=99,err=99)  iplevel             
          read(ioppar,*,end=99,err=99)  paplo               
          read(ioppar,*,end=99,err=99)  errplo              
          read(ioppar,*,end=99,err=99)  paxis               
          read(ioppar,*,end=99,err=99)  taxis               
          read(ioppar,*,end=99,err=99)  fitplo              
          read(ioppar,*,end=99,err=99)  log_x               
          read(ioppar,*,end=99,err=99)  log_y               
          read(ioppar,*,end=99,err=99)  log_x               
          read(ioppar,*,end=99,err=99)  log_y                  
          read(ioppar,*,end=99,err=99)  ibild               
          read(ioppar,*,end=99,err=99)  ifont               
          read(ioppar,*,end=99,err=99)  ax_text_scale       
          read(ioppar,*,end=99,err=99)  ax_tick_scale       
          read(ioppar,*,end=99,err=99)  thcline_thickness   
          read(ioppar,*,end=99,err=99)  datline_thickness   
          read(ioppar,*,end=99,err=99)  XLEG_DISTANCE       
          read(ioppar,*,end=99,err=99)  YLEG_DISTANCE       
          read(ioppar,*,end=99,err=99)  TLEG_DISTANCE_X     
          read(ioppar,*,end=99,err=99)  TLEG_DISTANCE_Y     
! to be modernized
          read(ioppar,*,end=99,err=99) nl
          read(ioppar,*,end=99,err=99) linetype(1:min(nl,size(linetype)))
          read(ioppar,*,end=99,err=99) nl
          read(ioppar,*,end=99,err=99) linewidth_scaling(1:min(nl,size(linewidth_scaling)))
          read(ioppar,*,end=99,err=99) nl
          read(ioppar,*,end=99,err=99) isymb(1:min(nl,size(isymb)))
          read(ioppar,*,end=99,err=99) nl
          read(ioppar,*,end=99,err=99) sysize_scaling(1:min(nl,size(sysize_scaling)))
          read(ioppar,*,end=99,err=99) nl
          read(ioppar,*,end=99,err=99) icolo(1:min(nl,size(icolo)))

99        continue
         close(ioppar)

       end subroutine readplotpar
   
       subroutine writeplotpar()
         implicit none
         integer :: ioppar
         
         open(newunit=ioppar, file=LAST_PPAR)
          write(ioppar,*)  xmin                ,'                      xmin    '    
          write(ioppar,*)  xmax                ,'                      xmax    ' 
          write(ioppar,*)  ymin                ,'                      ymin    ' 
          write(ioppar,*)  ymax                ,'                      ymax    ' 
          write(ioppar,*)  ptex                ,'                      text    ' 
          write(ioppar,*)  epscl               ,'                      epscl   ' 
          write(ioppar,*)  txsize              ,'                      txsize  ' 
          write(ioppar,*)  txsizt              ,'                      legsize ' 
          write(ioppar,*)  xtshft              ,'                      legx    ' 
          write(ioppar,*)  ytshft              ,'                      legy    ' 
          write(ioppar,*)  icol0               ,'                      color   ' 
          write(ioppar,*)  sysize              ,'                      sysize  ' 
          write(ioppar,*)  iplevel             ,'                      parlev  ' 
          write(ioppar,*)  paplo               ,'                      parplo  ' 
          write(ioppar,*)  errplo              ,'                      noerrors' 
          write(ioppar,*)  paxis               ,'                      axis    ' 
          write(ioppar,*)  taxis               ,'                      txaxis  ' 
          write(ioppar,*)  fitplo              ,'                      nofits  ' 
          write(ioppar,*)  log_x               ,'                      log_x   ' 
          write(ioppar,*)  log_y               ,'                      log_y   ' 
          write(ioppar,*)  log_x               ,'                      lin_x   ' 
          write(ioppar,*)  log_y               ,'                      lin_y   '    
          write(ioppar,*)  ibild               ,'                      #        '
          write(ioppar,*)  ifont               ,'                      font     '
          write(ioppar,*)  ax_text_scale       ,'                      axtxsize '
          write(ioppar,*)  ax_tick_scale       ,'                      axticlen '
          write(ioppar,*)  thcline_thickness   ,'                      flinewd  '
          write(ioppar,*)  datline_thickness   ,'                      dlinewd  '
          write(ioppar,*)  XLEG_DISTANCE       ,'                      xlegdist '
          write(ioppar,*)  YLEG_DISTANCE       ,'                      ylegdist '
          write(ioppar,*)  TLEG_DISTANCE_X     ,'                      tit_x    '
          write(ioppar,*)  TLEG_DISTANCE_Y     ,'                      tit_y    '
          write(ioppar,*) size(linetype)
          write(ioppar,'(20i5)') linetype
          write(ioppar,*) size(linewidth_scaling)
          write(ioppar,'(10f12.6)') linewidth_scaling
          write(ioppar,*) size(isymb)
          write(ioppar,'(20i5)') isymb
          write(ioppar,*) size(sysize_scaling)
          write(ioppar,'(10f12.6)') sysize_scaling
          write(ioppar,*) size(icolo)
          write(ioppar,'(20i5)') icolo
 
         close(ioppar)

       end subroutine writeplotpar
     
     
      END
!
!
       subroutine out_gli
!      ==================
!
       use new_com
       ! use cincoc
       use cdata
!       use outlev
       use selist
       implicit none

       character*12 infile
       logical*4    fileda
       integer i,j
       real xxx

! -- open the file -----------------------------------------------------
       if(inames.eq.0) then
         ierrs = 1
         write(6,*)'input(e1): filename is lacking!'
         return
       endif
!
       infile = vname(1)
       do i=1,8
         if(infile(i:i).eq.' ') infile(i:i)='_'
       enddo
       infile(9:12)='.dat'
       inquire(file=infile,exist=fileda)
!!     if(.not.fileda) then
!!       write(6,*)'input(e2): file ',infile,' does not exist!'
!!       ierrs = 3
!!       return
!!     endif

       if(nsel.gt.18 ) then
         write(6,*)'gli_out: too many selected items!'
         ierrs = 3
         return
       endif

       open(20,file=infile,status='UNKNOWN')
       write(20,'(A)')(name(isels(i)),i=1,nsel)
       do j=1,nwert(isels(1))
         write(20,'(18(1x,E13.6))')xwerte(j,isels(1)),                  &
     &        (ywerte(j,isels(i)),i=1,nsel)
         xxx = xwerte(j,isels(1))
         do i=1,nsel
           if( xxx-xwerte(j,isels(i)) .ne. 0.0 ) then
             write(6,*) 'Warning: x-value of ',isels(i),' differs!'
           endif
         enddo
       enddo

       close(20)
       write(6,*)'File: ',infile,' written!'
       return
      END
