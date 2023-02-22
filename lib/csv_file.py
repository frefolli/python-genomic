"""
    @does define CsvFile abstraction
"""

import logging
import os


class CsvFile:
    """
        @does define CsvFile abstraction
    """

    def __init__(self, path: str):
        """
            @does initialize CsvFile
        """
        self.path = path
        self.header = []
        self.body = []

    def __read_if_exists(self):
        """
            @does read content of csv file if exists in filesystem
        """
        if os.path.exists(self.path):
            self.__read()

    def __read(self):
        """
            @does read content of csvfile
        """
        with open(self.path, "r", encoding="utf-8") as file:
            self.__parse_read_lines(file.read().split("\n"))

    def __parse_read_lines(self, lines):
        """
            @does does parse content of lines to python data
        """
        if len(lines) > 0:
            self.header = lines[0].split(",")
            for line in lines[1:]:
                splitted_line = line.split(",")
                assert len(splitted_line) == len(self.header)
                self.body.append(splitted_line)

    def change_path(self, new_path: str):
        """
            @does change path of csv file
        """
        self.path = new_path

    def save(self):
        """
            @does save content of csv file to filesystem
        """
        assert len(self.header) > 0
        header_line = ",".join(self.header)
        body_lines = self.__dump_body_lines()
        self.__write("\n".join([header_line] + body_lines))
        logging.info("SAVED csvfile %s", self.path)

    def __dump_body_lines(self) -> list[str]:
        """
            @does dump python data to lines
        """
        body_lines = []
        for line in self.body:
            assert len(line) == len(self.header)
            body_lines.append(",".join([str(field) for field in line]))
        return body_lines

    def __write(self, lines: str):
        """
            @does write lines to filesystem
        """
        with open(self.path, "w", encoding="utf-8") as file:
            file.write(lines)

    def save_as(self, new_path: str):
        """
            @does rename the file and then save it
        """
        self.change_path(new_path)
        self.save()

    def rename_column(self, target_column: str, new_column_name: str):
        """
            @does rename a column
        """
        self.header[self.header.index(target_column)] = new_column_name

    def remove_column(self, target_column: str):
        """
            @does erase a column
        """
        index = self.header.index(target_column)
        self.header = self.header[0:index] + self.header[index+1:]
        self.body = [line[0:index] + line[index+1:] for line in self.body]

    def add_column(self, new_column_name: str):
        """
            @does add acolumn
        """
        assert new_column_name not in self.header
        self.header.append(new_column_name)
        self.body = [line + [""] for line in self.body]

    def rebase(self, column_list: list[str]):
        """
            @does reset whole content and set new header
        """
        assert len(column_list) > 0
        self.header = []
        self.body = []
        for column in column_list:
            self.add_column(column)
        logging.info("REBASED csvfile %s", self.path)

    def add_line(self, field_list: list):
        """
            @does add a line
        """
        assert len(field_list) == len(self.header)
        self.body.append(field_list)

    def remove_line(self, line_index: int):
        """
            @does remove a line
        """
        assert line_index >= 0
        assert line_index < len(self)
        self.body = self.body[0:line_index] + self.body[line_index+1:]

    def get_line(self, line_index: int):
        """
            @returns a line
        """
        assert line_index >= 0
        assert line_index < len(self)
        return self.body[line_index]

    def __len__(self) -> int:
        """
            @returns len
        """
        return len(self.body)

    def __enter__(self):
        """
            @does enter for with-as
        """
        self.__read_if_exists()
        logging.info("OPENED csvfile %s", self.path)
        return self

    def __exit__(self, exception_type,
                 exception_value, exception_traceback):
        """
            @does exit for with-as
        """
        logging.info("CLOSED csvfile %s", self.path)

    def __str__(self) -> str:
        """
            @returns string representation of CsvFile
        """
        return f"CsvFile(path = \"{self.path}\")"
