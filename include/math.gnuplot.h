#ifndef MATH_GNUPLOT_H
#define MATH_GNUPLOT_H

#include <cstddef>
#include <string>
#include <vector>

#include "defs.h"
#include "math.defs.h"

namespace gnu
{
    const char * const 	Snow            = "\'#FFFAFA\'";
    const char * const 	GhostWhite      = "\'#F8F8FF\'";
    const char * const 	WhiteSmoke      = "\'#F5F5F5\'";
    const char * const 	Gainsboro       = "\'#DCDCDC\'";
    const char * const 	FloralWhite     = "\'#FFFAF0\'";
    const char * const 	OldLace         = "\'#FDF5E6\'";
    const char * const 	Linen           = "\'#FAF0E6\'";
    const char * const 	AntiqueWhite	= "\'#FAEBD7\'";
    const char * const 	PapayaWhip      = "\'#FFEFD5\'";
    const char * const 	BlanchedAlmond	= "\'#FFEBCD\'";
    const char * const 	Bisque          = "\'#FFE4C4\'";
    const char * const 	PeachPuff       = "\'#FFDAB9\'";
    const char * const 	NavajoWhite     = "\'#FFDEAD\'";
    const char * const 	Moccasin        = "\'#FFE4B5\'";
    const char * const 	Cornsilk        = "\'#FFF8DC\'";
    const char * const 	Ivory           = "\'#FFFFF0\'";
    const char * const 	LemonChiffon    = "\'#FFFACD\'";
    const char * const 	Seashell        = "\'#FFF5EE\'";
    const char * const 	Honeydew        = "\'#F0FFF0\'";
    const char * const 	MintCream       = "\'#F5FFFA\'";
    const char * const 	Azure           = "\'#F0FFFF\'";
    const char * const 	AliceBlue       = "\'#F0F8FF\'";
    const char * const 	lavender        = "\'#E6E6FA\'";
    const char * const 	LavenderBlush   = "\'#FFF0F5\'";
    const char * const 	MistyRose       = "\'#FFE4E1\'";
    const char * const 	White           = "\'#FFFFFF\'";
    const char * const 	Black           = "\'#000000\'";
    const char * const 	DarkSlateGray   = "\'#2F4F4F\'";
    const char * const 	DimGrey         = "\'#696969\'";
    const char * const 	SlateGrey       = "\'#708090\'";
    const char * const 	LightSlateGray  = "\'#778899\'";
    const char * const 	Grey            = "\'#BEBEBE\'";
    const char * const 	LightGray       = "\'#D3D3D3\'";
    const char * const 	MidnightBlue    = "\'#191970\'";
    const char * const 	NavyBlue        = "\'#000080\'";
    const char * const 	CornflowerBlue  = "\'#6495ED\'";
    const char * const 	DarkSlateBlue   = "\'#483D8B\'";
    const char * const 	SlateBlue       = "\'#6A5ACD\'";
    const char * const 	MediumSlateBlue	= "\'#7B68EE\'";
    const char * const 	LightSlateBlue	= "\'#8470FF\'";
    const char * const 	MediumBlue      = "\'#0000CD\'";
    const char * const 	RoyalBlue       = "\'#4169E1\'";
    const char * const 	Blue            = "\'#0000FF\'";
    const char * const 	DodgerBlue      = "\'#1E90FF\'";
    const char * const 	DeepSkyBlue     = "\'#00BFFF\'";
    const char * const 	SkyBlue         = "\'#87CEEB\'";
    const char * const 	LightSkyBlue    = "\'#87CEFA\'";
    const char * const 	SteelBlue       = "\'#4682B4\'";
    const char * const 	LightSteelBlue  = "\'#B0C4DE\'";
    const char * const 	LightBlue       = "\'#ADD8E6\'";
    const char * const 	PowderBlue      = "\'#B0E0E6\'";
    const char * const 	PaleTurquoise   = "\'#AFEEEE\'";
    const char * const 	DarkTurquoise   = "\'#00CED1\'";
    const char * const 	MediumTurquoise = "\'#48D1CC\'";
    const char * const 	Turquoise       = "\'#40E0D0\'";
    const char * const 	Cyan            = "\'#00FFFF\'";
    const char * const 	LightCyan       = "\'#E0FFFF\'";
    const char * const 	CadetBlue       = "\'#5F9EA0\'";
    const char * const 	MediumAquamarine= "\'#66CDAA\'";
    const char * const 	Aquamarine      = "\'#7FFFD4\'";
    const char * const 	DarkGreen       = "\'#006400\'";
    const char * const 	DarkOliveGreen  = "\'#556B2F\'";
    const char * const 	DarkSeaGreen    = "\'#8FBC8F\'";
    const char * const 	SeaGreen        = "\'#2E8B57\'";
    const char * const 	MediumSeaGreen  = "\'#3CB371\'";
    const char * const 	LightSeaGreen   = "\'#20B2AA\'";
    const char * const 	PaleGreen       = "\'#98FB98\'";
    const char * const 	SpringGreen     = "\'#00FF7F\'";
    const char * const 	LawnGreen       = "\'#7CFC00\'";
    const char * const 	Green           = "\'#00FF00\'";
    const char * const 	Chartreuse      = "\'#7FFF00\'";
    const char * const 	MedSpringGreen  = "\'#00FA9A\'";
    const char * const 	GreenYellow     = "\'#ADFF2F\'";
    const char * const 	LimeGreen       = "\'#32CD32\'";
    const char * const 	YellowGreen     = "\'#9ACD32\'";
    const char * const 	ForestGreen     = "\'#228B22\'";
    const char * const 	OliveDrab       = "\'#6B8E23\'";
    const char * const 	DarkKhaki       = "\'#BDB76B\'";
    const char * const 	PaleGoldenrod   = "\'#EEE8AA\'";
    const char * const 	LtGoldenrodYello= "\'#FAFAD2\'";
    const char * const 	LightYellow     = "\'#FFFFE0\'";
    const char * const 	Yellow          = "\'#FFFF00\'";
    const char * const 	Gold            = "\'#FFD700\'";
    const char * const 	LightGoldenrod  = "\'#EEDD82\'";
    const char * const 	goldenrod       = "\'#DAA520\'";
    const char * const 	DarkGoldenrod   = "\'#B8860B\'";
    const char * const 	RosyBrown       = "\'#BC8F8F\'";
    const char * const 	IndianRed       = "\'#CD5C5C\'";
    const char * const 	SaddleBrown     = "\'#8B4513\'";
    const char * const 	Sienna          = "\'#A0522D\'";
    const char * const 	Peru            = "\'#CD853F\'";
    const char * const 	Burlywood       = "\'#DEB887\'";
    const char * const 	Beige           = "\'#F5F5DC\'";
    const char * const 	Wheat           = "\'#F5DEB3\'";
    const char * const 	SandyBrown      = "\'#F4A460\'";
    const char * const 	Tan             = "\'#D2B48C\'";
    const char * const 	Chocolate       = "\'#D2691E\'";
    const char * const 	Firebrick       = "\'#B22222\'";
    const char * const 	Brown           = "\'#A52A2A\'";
    const char * const 	DarkSalmon      = "\'#E9967A\'";
    const char * const 	Salmon          = "\'#FA8072\'";
    const char * const 	LightSalmon     = "\'#FFA07A\'";
    const char * const 	Orange          = "\'#FFA500\'";
    const char * const 	DarkOrange      = "\'#FF8C00\'";
    const char * const 	Coral           = "\'#FF7F50\'";
    const char * const 	LightCoral      = "\'#F08080\'";
    const char * const 	Tomato          = "\'#FF6347\'";
    const char * const 	OrangeRed       = "\'#FF4500\'";
    const char * const 	Red             = "\'#FF0000\'";
    const char * const 	HotPink         = "\'#FF69B4\'";
    const char * const 	DeepPink        = "\'#FF1493\'";
    const char * const 	Pink            = "\'#FFC0CB\'";
    const char * const 	LightPink       = "\'#FFB6C1\'";
    const char * const 	PaleVioletRed   = "\'#DB7093\'";
    const char * const 	Maroon          = "\'#B03060\'";
    const char * const 	MediumVioletRed = "\'#C71585\'";
    const char * const 	VioletRed       = "\'#D02090\'";
    const char * const 	Magenta         = "\'#FF00FF\'";
    const char * const 	Violet          = "\'#EE82EE\'";
    const char * const 	Plum            = "\'#DDA0DD\'";
    const char * const 	Orchid          = "\'#DA70D6\'";
    const char * const 	MediumOrchid    = "\'#BA55D3\'";
    const char * const 	DarkOrchid      = "\'#9932CC\'";
    const char * const 	DarkViolet      = "\'#9400D3\'";
    const char * const 	BlueViolet      = "\'#8A2BE2\'";
    const char * const 	Purple          = "\'#A020F0\'";
    const char * const 	MediumPurple    = "\'#9370DB\'";
    const char * const 	Thistle         = "\'#D8BFD8\'";
    const char * const 	Snow1           = "\'#FFFAFA\'";
    const char * const 	Snow2           = "\'#EEE9E9\'";
    const char * const 	Snow3           = "\'#CDC9C9\'";
    const char * const 	Snow4           = "\'#8B8989\'";
    const char * const 	Seashell1       = "\'#FFF5EE\'";
    const char * const 	Seashell2       = "\'#EEE5DE\'";
    const char * const 	Seashell3       = "\'#CDC5BF\'";
    const char * const 	Seashell4       = "\'#8B8682\'";
    const char * const 	AntiqueWhite1   = "\'#FFEFDB\'";
    const char * const 	AntiqueWhite2   = "\'#EEDFCC\'";
    const char * const 	AntiqueWhite3   = "\'#CDC0B0\'";
    const char * const 	AntiqueWhite4   = "\'#8B8378\'";
    const char * const 	Bisque1         = "\'#FFE4C4\'";
    const char * const 	Bisque2         = "\'#EED5B7\'";
    const char * const 	Bisque3         = "\'#CDB79E\'";
    const char * const 	Bisque4         = "\'#8B7D6B\'";
    const char * const 	PeachPuff1      = "\'#FFDAB9\'";
    const char * const 	PeachPuff2      = "\'#EECBAD\'";
    const char * const 	PeachPuff3      = "\'#CDAF95\'";
    const char * const 	PeachPuff4      = "\'#8B7765\'";
    const char * const 	NavajoWhite1    = "\'#FFDEAD\'";
    const char * const 	NavajoWhite2    = "\'#EECFA1\'";
    const char * const 	NavajoWhite3    = "\'#CDB38B\'";
    const char * const 	NavajoWhite4    = "\'#8B795E\'";
    const char * const 	LemonChiffon1   = "\'#FFFACD\'";
    const char * const 	LemonChiffon2   = "\'#EEE9BF\'";
    const char * const 	LemonChiffon3   = "\'#CDC9A5\'";
    const char * const 	LemonChiffon4	= "\'#8B8970\'";
    const char * const 	Cornsilk1	= "\'#FFF8DC\'";
    const char * const 	Cornsilk2	= "\'#EEE8CD\'";
    const char * const 	Cornsilk3	= "\'#CDC8B1\'";
    const char * const 	Cornsilk4	= "\'#8B8878\'";
    const char * const 	vory1	= "\'#FFFFF0\'";
    const char * const 	vory2	= "\'#EEEEE0\'";
    const char * const 	Ivory3	= "\'#CDCDC1\'";
    const char * const 	Ivory4	= "\'#8B8B83\'";
    const char * const 	Honeydew1	= "\'#F0FFF0\'";
    const char * const 	Honeydew2	= "\'#E0EEE0\'";
    const char * const 	Honeydew3	= "\'#C1CDC1\'";
    const char * const 	Honeydew4	= "\'#838B83\'";
    const char * const 	LavenderBlush1	= "\'#FFF0F5\'";
    const char * const 	LavenderBlush2	= "\'#EEE0E5\'";
    const char * const 	LavenderBlush3	= "\'#CDC1C5\'";
    const char * const 	LavenderBlush4	= "\'#8B8386\'";
    const char * const 	MistyRose1	= "\'#FFE4E1\'";
    const char * const 	MistyRose2	= "\'#EED5D2\'";
    const char * const 	MistyRose3	= "\'#CDB7B5\'";
    const char * const 	MistyRose4	= "\'#8B7D7B\'";
    const char * const 	Azure1	= "\'#F0FFFF\'";
    const char * const 	Azure2	= "\'#E0EEEE\'";
    const char * const 	Azure3	= "\'#C1CDCD\'";
    const char * const 	Azure4	= "\'#838B8B\'";
    const char * const 	SlateBlue1	= "\'#836FFF\'";
    const char * const 	SlateBlue2	= "\'#7A67EE\'";
    const char * const 	SlateBlue3	= "\'#6959CD\'";
    const char * const 	SlateBlue4	= "\'#473C8B\'";
    const char * const 	RoyalBlue1	= "\'#4876FF\'";
    const char * const 	RoyalBlue2	= "\'#436EEE\'";
    const char * const 	RoyalBlue3	= "\'#3A5FCD\'";
    const char * const 	RoyalBlue4	= "\'#27408B\'";
    const char * const 	Blue1	= "\'#0000FF\'";
    const char * const 	Blue2	= "\'#0000EE\'";
    const char * const 	Blue3	= "\'#0000CD\'";
    const char * const 	Blue4	= "\'#00008B\'";
    const char * const 	DodgerBlue1	= "\'#1E90FF\'";
    const char * const 	DodgerBlue2	= "\'#1C86EE\'";
    const char * const 	DodgerBlue3	= "\'#1874CD\'";
    const char * const 	DodgerBlue4	= "\'#104E8B\'";
    const char * const 	SteelBlue1	= "\'#63B8FF\'";
    const char * const 	SteelBlue2	= "\'#5CACEE\'";
    const char * const 	SteelBlue3	= "\'#4F94CD\'";
    const char * const 	SteelBlue4	= "\'#36648B\'";
    const char * const 	DeepSkyBlue1	= "\'#00BFFF\'";
    const char * const 	DeepSkyBlue2	= "\'#00B2EE\'";
    const char * const 	DeepSkyBlue3	= "\'#009ACD\'";
    const char * const 	DeepSkyBlue4	= "\'#00688B\'";
    const char * const 	SkyBlue1	= "\'#87CEFF\'";
    const char * const 	SkyBlue2	= "\'#7EC0EE\'";
    const char * const 	SkyBlue3	= "\'#6CA6CD\'";
    const char * const 	SkyBlue4	= "\'#4A708B\'";
    const char * const 	LightSkyBlue1	= "\'#B0E2FF\'";
    const char * const 	LightSkyBlue2	= "\'#A4D3EE\'";
    const char * const 	LightSkyBlue3	= "\'#8DB6CD\'";
    const char * const 	LightSkyBlue4	= "\'#607B8B\'";
    const char * const 	SlateGray1	= "\'#C6E2FF\'";
    const char * const 	SlateGray2	= "\'#B9D3EE\'";
    const char * const 	SlateGray3	= "\'#9FB6CD\'";
    const char * const 	SlateGray4	= "\'#6C7B8B\'";
    const char * const 	LightSteelBlue1	= "\'#CAE1FF\'";
    const char * const 	LightSteelBlue2	= "\'#BCD2EE\'";
    const char * const 	LightSteelBlue3	= "\'#A2B5CD\'";
    const char * const 	LightSteelBlue4	= "\'#6E7B8B\'";
    const char * const 	LightBlue1	= "\'#BFEFFF\'";
    const char * const 	LightBlue2	= "\'#B2DFEE\'";
    const char * const 	LightBlue3	= "\'#9AC0CD\'";
    const char * const 	LightBlue4	= "\'#68838B\'";
    const char * const 	LightCyan1	= "\'#E0FFFF\'";
    const char * const 	LightCyan2	= "\'#D1EEEE\'";
    const char * const 	LightCyan3	= "\'#B4CDCD\'";
    const char * const 	LightCyan4	= "\'#7A8B8B\'";
    const char * const 	PaleTurquoise1	= "\'#BBFFFF\'";
    const char * const 	PaleTurquoise2	= "\'#AEEEEE\'";
    const char * const 	PaleTurquoise3	= "\'#96CDCD\'";
    const char * const 	PaleTurquoise4	= "\'#668B8B\'";
    const char * const 	CadetBlue1	= "\'#98F5FF\'";
    const char * const 	CadetBlue2	= "\'#8EE5EE\'";
    const char * const 	CadetBlue3	= "\'#7AC5CD\'";
    const char * const 	CadetBlue4	= "\'#53868B\'";
    const char * const 	Turquoise1	= "\'#00F5FF\'";
    const char * const 	Turquoise2	= "\'#00E5EE\'";
    const char * const 	Turquoise3	= "\'#00C5CD\'";
    const char * const 	Turquoise4	= "\'#00868B\'";
    const char * const 	Cyan1	= "\'#00FFFF\'";
    const char * const 	Cyan2	= "\'#00EEEE\'";
    const char * const 	Cyan3	= "\'#00CDCD\'";
    const char * const 	Cyan4	= "\'#008B8B\'";
    const char * const 	DarkSlateGray1	= "\'#97FFFF\'";
    const char * const 	DarkSlateGray2	= "\'#8DEEEE\'";
    const char * const 	DarkSlateGray3	= "\'#79CDCD\'";
    const char * const 	DarkSlateGray4	= "\'#528B8B\'";
    const char * const 	Aquamarine1	= "\'#7FFFD4\'";
    const char * const 	Aquamarine2	= "\'#76EEC6\'";
    const char * const 	Aquamarine3	= "\'#66CDAA\'";
    const char * const 	Aquamarine4	= "\'#458B74\'";
    const char * const 	DarkSeaGreen1	= "\'#C1FFC1\'";
    const char * const 	DarkSeaGreen2	= "\'#B4EEB4\'";
    const char * const 	DarkSeaGreen3	= "\'#9BCD9B\'";
    const char * const 	DarkSeaGreen4	= "\'#698B69\'";
    const char * const 	SeaGreen1	= "\'#54FF9F\'";
    const char * const 	SeaGreen2	= "\'#4EEE94\'";
    const char * const 	SeaGreen3	= "\'#43CD80\'";
    const char * const 	SeaGreen4	= "\'#2E8B57\'";
    const char * const 	PaleGreen1	= "\'#9AFF9A\'";
    const char * const 	PaleGreen2	= "\'#90EE90\'";
    const char * const 	PaleGreen3	= "\'#7CCD7C\'";
    const char * const 	PaleGreen4	= "\'#548B54\'";
    const char * const 	SpringGreen1	= "\'#00FF7F\'";
    const char * const 	SpringGreen2	= "\'#00EE76\'";
    const char * const 	SpringGreen3	= "\'#00CD66\'";
    const char * const 	SpringGreen4	= "\'#008B45\'";
    const char * const 	Green1	= "\'#00FF00\'";
    const char * const 	Green2	= "\'#00EE00\'";
    const char * const 	Green3	= "\'#00CD00\'";
    const char * const 	Green4	= "\'#008B00\'";
    const char * const 	Chartreuse1	= "\'#7FFF00\'";
    const char * const 	Chartreuse2	= "\'#76EE00\'";
    const char * const 	Chartreuse3	= "\'#66CD00\'";
    const char * const 	Chartreuse4	= "\'#458B00\'";
    const char * const 	OliveDrab1	= "\'#C0FF3E\'";
    const char * const 	OliveDrab2	= "\'#B3EE3A\'";
    const char * const 	OliveDrab3	= "\'#9ACD32\'";
    const char * const 	OliveDrab4	= "\'#698B22\'";
    const char * const 	DarkOliveGreen1	= "\'#CAFF70\'";
    const char * const 	DarkOliveGreen2	= "\'#BCEE68\'";
    const char * const 	DarkOliveGreen3	= "\'#A2CD5A\'";
    const char * const 	DarkOliveGreen4	= "\'#6E8B3D\'";
    const char * const 	Khaki1	= "\'#FFF68F\'";
    const char * const 	Khaki2	= "\'#EEE685\'";
    const char * const 	Khaki3	= "\'#CDC673\'";
    const char * const 	Khaki4	= "\'#8B864E\'";
    const char * const 	LightGoldenrod1	= "\'#FFEC8B\'";
    const char * const 	LightGoldenrod2	= "\'#EEDC82\'";
    const char * const 	LightGoldenrod3	= "\'#CDBE70\'";
    const char * const 	LightGoldenrod4	= "\'#8B814C\'";
    const char * const 	LightYellow1	= "\'#FFFFE0\'";
    const char * const 	LightYellow2	= "\'#EEEED1\'";
    const char * const 	LightYellow3	= "\'#CDCDB4\'";
    const char * const 	LightYellow4	= "\'#8B8B7A\'";
    const char * const 	Yellow1	= "\'#FFFF00\'";
    const char * const 	Yellow2	= "\'#EEEE00\'";
    const char * const 	Yellow3	= "\'#CDCD00\'";
    const char * const 	Yellow4	= "\'#8B8B00\'";
    const char * const 	Gold1	= "\'#FFD700\'";
    const char * const 	Gold2	= "\'#EEC900\'";
    const char * const 	Gold3	= "\'#CDAD00\'";
    const char * const 	Gold4	= "\'#8B7500\'";
    const char * const 	Goldenrod1	= "\'#FFC125\'";
    const char * const 	Goldenrod2	= "\'#EEB422\'";
    const char * const 	Goldenrod3	= "\'#CD9B1D\'";
    const char * const 	Goldenrod4	= "\'#8B6914\'";
    const char * const 	DarkGoldenrod1	= "\'#FFB90F\'";
    const char * const 	DarkGoldenrod2	= "\'#EEAD0E\'";
    const char * const 	DarkGoldenrod3	= "\'#CD950C\'";
    const char * const 	DarkGoldenrod4	= "\'#8B658B\'";
    const char * const 	RosyBrown1	= "\'#FFC1C1\'";
    const char * const 	RosyBrown2	= "\'#EEB4B4\'";
    const char * const 	RosyBrown3	= "\'#CD9B9B\'";
    const char * const 	RosyBrown4	= "\'#8B6969\'";
    const char * const 	IndianRed1	= "\'#FF6A6A\'";
    const char * const 	IndianRed2	= "\'#EE6363\'";
    const char * const 	IndianRed3	= "\'#CD5555\'";
    const char * const 	IndianRed4	= "\'#8B3A3A\'";
    const char * const 	Sienna1	= "\'#FF8247\'";
    const char * const 	Sienna2	= "\'#EE7942\'";
    const char * const 	Sienna3	= "\'#CD6839\'";
    const char * const 	Sienna4	= "\'#8B4726\'";
    const char * const 	Burlywood1	= "\'#FFD39B\'";
    const char * const 	Burlywood2	= "\'#EEC591\'";
    const char * const 	Burlywood3	= "\'#CDAA7D\'";
    const char * const 	Burlywood4	= "\'#8B7355\'";
    const char * const 	Wheat1	= "\'#FFE7BA\'";
    const char * const 	Wheat2	= "\'#EED8AE\'";
    const char * const 	Wheat3	= "\'#CDBA96\'";
    const char * const 	Wheat4	= "\'#8B7E66\'";
    const char * const 	Tan1	= "\'#FFA54F\'";
    const char * const 	Tan2	= "\'#EE9A49\'";
    const char * const 	Tan3	= "\'#CD853F\'";
    const char * const 	Tan4	= "\'#8B5A2B\'";
    const char * const 	Chocolate1	= "\'#FF7F24\'";
    const char * const 	Chocolate2	= "\'#EE7621\'";
    const char * const 	Chocolate3	= "\'#CD661D\'";
    const char * const 	Chocolate4	= "\'#8B4513\'";
    const char * const 	Firebrick1	= "\'#FF3030\'";
    const char * const 	Firebrick2	= "\'#EE2C2C\'";
    const char * const 	Firebrick3	= "\'#CD2626\'";
    const char * const 	Firebrick4	= "\'#8B1A1A\'";
    const char * const 	Brown1	= "\'#FF4040\'";
    const char * const 	Brown2	= "\'#EE3B3B\'";
    const char * const 	Brown3	= "\'#CD3333\'";
    const char * const 	Brown4	= "\'#8B2323\'";
    const char * const 	Salmon1	= "\'#FF8C69\'";
    const char * const 	Salmon2	= "\'#EE8262\'";
    const char * const 	Salmon3	= "\'#CD7054\'";
    const char * const 	Salmon4	= "\'#8B4C39\'";
    const char * const 	LightSalmon1	= "\'#FFA07A\'";
    const char * const 	LightSalmon2	= "\'#EE9572\'";
    const char * const 	LightSalmon3	= "\'#CD8162\'";
    const char * const 	LightSalmon4	= "\'#8B5742\'";
    const char * const 	Orange1	= "\'#FFA500\'";
    const char * const 	Orange2	= "\'#EE9A00\'";
    const char * const 	Orange3	= "\'#CD8500\'";
    const char * const 	Orange4	= "\'#8B5A00\'";
    const char * const 	DarkOrange1	= "\'#FF7F00\'";
    const char * const 	DarkOrange2	= "\'#EE7600\'";
    const char * const 	DarkOrange3	= "\'#CD6600\'";
    const char * const 	DarkOrange4	= "\'#8B4500\'";
    const char * const 	Coral1	= "\'#FF7256\'";
    const char * const 	Coral2	= "\'#EE6A50\'";
    const char * const 	Coral3	= "\'#CD5B45\'";
    const char * const 	Coral4	= "\'#8B3E2F\'";
    const char * const 	Tomato1	= "\'#FF6347\'";
    const char * const 	Tomato2	= "\'#EE5C42\'";
    const char * const 	Tomato3	= "\'#CD4F39\'";
    const char * const 	Tomato4	= "\'#8B3626\'";
    const char * const 	OrangeRed1	= "\'#FF4500\'";
    const char * const 	OrangeRed2	= "\'#EE4000\'";
    const char * const 	OrangeRed3	= "\'#CD3700\'";
    const char * const 	OrangeRed4	= "\'#8B2500\'";
    const char * const 	Red1	= "\'#FF0000\'";
    const char * const 	Red2	= "\'#EE0000\'";
    const char * const 	Red3	= "\'#CD0000\'";
    const char * const 	Red4	= "\'#8B0000\'";
    const char * const 	DeepPink1	= "\'#FF1493\'";
    const char * const 	DeepPink2	= "\'#EE1289\'";
    const char * const 	DeepPink3	= "\'#CD1076\'";
    const char * const 	DeepPink4	= "\'#8B0A50\'";
    const char * const 	HotPink1	= "\'#FF6EB4\'";
    const char * const 	HotPink2	= "\'#EE6AA7\'";
    const char * const 	HotPink3	= "\'#CD6090\'";
    const char * const 	HotPink4	= "\'#8B3A62\'";
    const char * const 	Pink1	= "\'#FFB5C5\'";
    const char * const 	Pink2	= "\'#EEA9B8\'";
    const char * const 	Pink3	= "\'#CD919E\'";
    const char * const 	Pink4	= "\'#8B636C\'";
    const char * const 	LightPink1	= "\'#FFAEB9\'";
    const char * const 	LightPink2	= "\'#EEA2AD\'";
    const char * const 	LightPink3	= "\'#CD8C95\'";
    const char * const 	LightPink4	= "\'#8B5F65\'";
    const char * const 	PaleVioletRed1	= "\'#FF82AB\'";
    const char * const 	PaleVioletRed2	= "\'#EE799F\'";
    const char * const 	PaleVioletRed3	= "\'#CD6889\'";
    const char * const 	PaleVioletRed4	= "\'#8B475D\'";
    const char * const 	Maroon1	= "\'#FF34B3\'";
    const char * const 	Maroon2	= "\'#EE30A7\'";
    const char * const 	Maroon3	= "\'#CD2990\'";
    const char * const 	Maroon4	= "\'#8B1C62\'";
    const char * const 	VioletRed1	= "\'#FF3E96\'";
    const char * const 	VioletRed2	= "\'#EE3A8C\'";
    const char * const 	VioletRed3	= "\'#CD3278\'";
    const char * const 	VioletRed4	= "\'#8B2252\'";
    const char * const 	Magenta1	= "\'#FF00FF\'";
    const char * const 	Magenta2	= "\'#EE00EE\'";
    const char * const 	Magenta3	= "\'#CD00CD\'";
    const char * const 	Magenta4	= "\'#8B008B\'";
    const char * const 	Orchid1	= "\'#FF83FA\'";
    const char * const 	Orchid2	= "\'#EE7AE9\'";
    const char * const 	Orchid3	= "\'#CD69C9\'";
    const char * const 	Orchid4	= "\'#8B4789\'";
    const char * const 	Plum1	= "\'#FFBBFF\'";
    const char * const 	Plum2	= "\'#EEAEEE\'";
    const char * const 	Plum3	= "\'#CD96CD\'";
    const char * const 	Plum4	= "\'#8B668B\'";
    const char * const 	MediumOrchid1	= "\'#E066FF\'";
    const char * const 	MediumOrchid2	= "\'#D15FEE\'";
    const char * const 	MediumOrchid3	= "\'#B452CD\'";
    const char * const 	MediumOrchid4	= "\'#7A378B\'";
    const char * const 	DarkOrchid1	= "\'#BF3EFF\'";
    const char * const 	DarkOrchid2	= "\'#B23AEE\'";
    const char * const 	DarkOrchid3	= "\'#9A32CD\'";
    const char * const 	DarkOrchid4	= "\'#68228B\'";
    const char * const 	Purple1	= "\'#9B30FF\'";
    const char * const 	Purple2	= "\'#912CEE\'";
    const char * const 	Purple3	= "\'#7D26CD\'";
    const char * const 	Purple4	= "\'#551A8B\'";
    const char * const 	MediumPurple1	= "\'#AB82FF\'";
    const char * const 	MediumPurple2	= "\'#9F79EE\'";
    const char * const 	MediumPurple3	= "\'#8968CD\'";
    const char * const 	MediumPurple4	= "\'#5D478B\'";
    const char * const 	Thistle1	= "\'#FFE1FF\'";
    const char * const 	Thistle2	= "\'#EED2EE\'";
    const char * const 	Thistle3	= "\'#CDB5CD\'";
    const char * const 	Thistle4	= "\'#8B7B8B\'";
    const char * const 	grey11	= "\'#1C1C1C\'";
    const char * const 	grey21	= "\'#363636\'";
    const char * const 	grey31	= "\'#4F4F4F\'";
    const char * const 	grey41	= "\'#696969\'";
    const char * const 	grey51	= "\'#828282\'";
    const char * const 	grey61	= "\'#9C9C9C\'";
    const char * const 	grey71	= "\'#B5B5B5\'";
    const char * const 	gray81	= "\'#CFCFCF\'";
    const char * const 	gray91	= "\'#E8E8E8\'";
    const char * const 	DarkGrey	= "\'#A9A9A9\'";
    const char * const 	DarkBlue	= "\'#00008B\'";
    const char * const 	DarkCyan	= "\'#008B8B\'";
    const char * const 	DarkMagenta	= "\'#8B008B\'";
    const char * const 	DarkRed	= "\'#8B0000\'";
    const char * const 	LightGreen	= "\'#90EE90\'";
} // namespace gnu

struct gnuplot_packet : public std::string
{
    gnuplot_packet() = default;
    gnuplot_packet(const char *s) { this->append(s); }
    gnuplot_packet(const std::string &s) { this->append(s); }
    gnuplot_packet(const gnuplot_packet &) = default;

    bool empty() const { return std::string::empty(); }

    gnuplot_packet &operator<< (const char *s) { this->append(s); return *this; }
    gnuplot_packet &operator<< (const std::string &s) { this->append(s); return *this; }
    gnuplot_packet &operator<< (char s) { this->append(1, s); return *this; }
    gnuplot_packet &operator<< (const gnuplot_packet &p) { this->append(p); return *this; }

    gnuplot_packet operator()(const char *s) { gnuplot_packet g(*this); g.append(s); return g; }
    gnuplot_packet operator()(const std::string &s) { gnuplot_packet g(*this); g.append(s); return g; }
    gnuplot_packet operator()(char s) { gnuplot_packet g(*this); g.append(1, s); return g; }

    template <typename T> gnuplot_packet &operator<< (T s) { this->append(std::to_string(s)); return *this; }
    template <typename T> gnuplot_packet &operator()(T s) { this->append(std::to_string(s)); return *this; }

    gnuplot_packet &wcircles(double r=0.02, const char *c="black", double t=1.)
    {
        this->append(":(").append(std::to_string(r)).append(") ")
            .append(" w circles fillcolor ").append(1,'\'').append(c).append(1,'\'')
            .append(" fs noborder transparent solid ").append(std::to_string(t));
        return *this;
    }
    gnuplot_packet &wlines(unsigned w=1, const char *c="black", int t=1)
    {
        this->append(" w line ").append(" lw ").append(std::to_string(w))
            .append(" lc ").append(1,'\'').append(c).append(1,'\'')
            .append(" lt ").append(std::to_string(t));
        return *this;
    }
};

template <typename T>
inline gnuplot_packet operator+(const gnuplot_packet &a, const T &b)
    { return gnuplot_packet(a) << b; }

template <typename T>
inline gnuplot_packet operator+(const T &a, const gnuplot_packet &b)
    { return gnuplot_packet(a) << b; }

class gnuplot
{
    FILE *pipe = nullptr;

public:

    gnuplot() { pipe = popen ("gnuplot -p", "w"); }
    ~gnuplot() { fflush(pipe); pclose(pipe); }
    gnuplot &operator<< (const gnuplot_packet &s) { fprintf(pipe, "%s\n", s.c_str()); return *this; }
};

/// стандартное имя временного файла
#define TEMP_FILE "pattern.temp"

template <unsigned n>
/*!
 * \brief подготовка файлов данных для вызова gnuplot
 * \note все файлы имеют стандартный префикс имени TEMP_FILE
 * \param patterns данные
 */
inline void gnuplot_prepare(const std::vector<curve<n>> &patterns)
{
    size_t np = patterns.size();

    for (size_t ipattern=0; ipattern<np; ipattern++)
    {
        FILE * temp = fopen((TEMP_FILE + std::to_string(ipattern)).c_str(), "w");
        for (auto e : patterns[ipattern])
        {
            for (unsigned i=0; i<n; i++) fprintf(temp, "%4.3lf ", e[i]);
            fprintf(temp, "\n");
        }
        fclose(temp);
    }
}

template <unsigned n>
/*!
 * \brief вставить данные в текущий график
 * \param pipe сценарий gnuplot
 * \param patterns данные по кривым
 * \param options опции отрисоки кривых
 * \param istart номер первого графика
 * \param iend номер за последним включаемым графиком
 */
inline void gnuplot_insert_plot(gnuplot &pipe, const std::vector<curve<n>> &patterns,
    const std::vector<gnuplot_packet> &options, Unsigned istart, Unsigned iend)
{
    assert(patterns.size() == options.size());
    assert(iend <= patterns.size());
    if (istart == iend) return;

    gnuplot_packet gnu_plot("plot");
    for (size_t ipattern=istart; ipattern<iend; ipattern++)
    {
        if (ipattern != istart) gnu_plot << ",";
        std::string tmpfile = std::string("'") + TEMP_FILE + std::to_string(ipattern) + "' ";
        gnu_plot << tmpfile << options[ipattern];
    }
    pipe << gnu_plot;
}

#endif // MATH_GNUPLOT_H