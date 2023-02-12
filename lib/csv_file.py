import os

class CsvFile:
    def __init__(self, path : str):
        self.path = path
        self.header = []
        self.body = []

    def __read_if_exists(self):
        if os.path.exists(self.path):
            self.__read()
    
    def __read(self):
        with open(self.path, "r") as file:
            self.__parse_read_lines(file.read().split("\n"))
    
    def __parse_read_lines(self, lines):
        if len(lines) > 0:
            self.header = lines[0].split(",")
            for line in lines[1:]:
                splitted_line = line.split(",")
                assert len(splitted_line) == len(self.header)
                self.body.append(splitted_line)
    
    def change_path(self, new_path : str):
        self.path = new_path
    
    def save(self):
        assert len(self.header) > 0
        header_line = ",".join(self.header)
        body_lines = self.__dump_body_lines()
        self.__write("\n".join([header_line] + body_lines))
    
    def __dump_body_lines(self) -> list[str]:
        body_lines = []
        for line in self.body:
            assert len(line) == len(self.header)
            body_lines.append(",".join([str(field) for field in line]))
        return body_lines
    
    def __write(self, lines : str):
        with open(self.path, "w") as file:
            file.write(lines)

    def save_as(self, new_path : str):
        self.change_path(new_path)
        self.save()
    
    def rename_column(self, target_column : str, new_column_name : str):
        self.header[self.header.index(target_column)] = new_column_name
    
    def remove_column(self, target_column : str):
        index = self.header.index(target_column)
        self.header = self.header[0:index] + self.header[index+1:]
        self.body = [line[0:index] + line[index+1:] for line in self.body]

    def add_column(self, new_column_name : str):
        assert new_column_name not in self.header
        self.header.append(new_column_name)
        self.body = [line + [""] for line in self.body]
    
    def rebase(self, column_list : list[str]):
        assert len(column_list) > 0
        self.header = []
        self.body = []
        for column in column_list: self.add_column(column)
    
    def add_line(self, field_list : list):
        assert len(field_list) == len(self.header)
        self.body.append(field_list)
    
    def remove_line(self, line_index : int):
        assert line_index >= 0
        assert line_index < self.__len__()
        self.body = self.body[0:line_index] + self.body[line_index+1:]
    
    def get_line(self, line_index : int):
        assert line_index >= 0
        assert line_index < self.__len__()
        return self.body[line_index]
    
    def __len__(self):
        return len(self.body)

    def __enter__(self):
        self.__read_if_exists()
        return self
    
    def __exit__(self, type, value, tb):
        # every change to referred file is explicit
        pass

    def __str__(self):
        return f"CsvFile(path = \"{self.path}\")"