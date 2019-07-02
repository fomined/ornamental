#ifndef TEST_INPUT_FILE_H
#define TEST_INPUT_FILE_H

#include <cstddef>
#include <string>
#include <map>

/*!
 * \brief тестирование файла на реальные и потенциальные ошибки
 * \param filename имя файла
 */
void test_file(const std::string &filename);

/*!
 * \brief тестирование всех файлов в заданной директории на реальные и потенциальные ошибки
 * \param filename имя директории
 */
void test_directory(const std::string &dirname);

#endif // TEST_INPUT_FILE_H
