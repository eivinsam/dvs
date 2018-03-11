#include <dvs.h>

int main()
{
	dvs::Space<float> R3(3);

	const auto foo = R3(1, 2, 3);
	const auto bar = foo.space()(4, 5, 6);

	const auto baz = R3(0, 1, 2) - R3(2, 1, 0);

	return 0;
}
