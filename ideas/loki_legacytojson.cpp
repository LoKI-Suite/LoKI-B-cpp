#include "ideas/LegacyToJSON.h"

int main(int argc, const char* argv[])
{
	loki::json_type json;
	if (argc==2)
	{
		json = loki::legacyToJSON(argv[1]);
	}
	else
	{
		json = loki::legacyToJSON(std::cin);
	}
	std::cout << json.dump(2) << std::endl;
	return 0;
}
