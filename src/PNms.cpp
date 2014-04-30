#include "PNms.h"

void PNms::readPNms(FileNode pNmsNode)
{
	type = (string)pNmsNode["type"];
	overlap = pNmsNode["overlap"];
	ovrDnm = (string)pNmsNode["ovrDnm"];
}
