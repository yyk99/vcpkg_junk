#include "unity.h"
#include "calculator.h"

void setUp(void) {
    // This is run before each test
}

void tearDown(void) {
    // This is run after each test
}

void test_add_positive_numbers(void) {
    TEST_ASSERT_EQUAL(5, add(2, 3));
    TEST_ASSERT_EQUAL(10, add(7, 3));
}

void test_add_negative_numbers(void) {
    TEST_ASSERT_EQUAL(-5, add(-2, -3));
    TEST_ASSERT_EQUAL(0, add(-5, 5));
}

void test_subtract(void) {
    TEST_ASSERT_EQUAL(2, subtract(5, 3));
    TEST_ASSERT_EQUAL(-2, subtract(3, 5));
    TEST_ASSERT_EQUAL(0, subtract(5, 5));
}

void test_multiply(void) {
    TEST_ASSERT_EQUAL(6, multiply(2, 3));
    TEST_ASSERT_EQUAL(-10, multiply(-2, 5));
    TEST_ASSERT_EQUAL(0, multiply(0, 5));
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_add_positive_numbers);
    RUN_TEST(test_add_negative_numbers);
    RUN_TEST(test_subtract);
    RUN_TEST(test_multiply);

    return UNITY_END();
}
