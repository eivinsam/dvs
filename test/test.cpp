#include <dvs.h>

int main()
{
	dvs::Space<float> R3(3);

	const auto foo = R3(1, 2, 3);
	const auto bar = foo.space()(4, 5, 6);


	return 0;
}
