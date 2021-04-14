//
// sql.h
//
// useful sql functions
//
// Copyright Peter Andrews 2021 @ CSHL
//

#ifndef UTILITY_SQL_H_
#define UTILITY_SQL_H_

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wdelete-non-abstract-non-virtual-dtor"
#include <mysqlx/xdevapi.h>
#pragma GCC diagnostic pop

#include <algorithm>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "html.h"
#include "paastrings.h"
#include "private.h"
#include "utility.h"

namespace paa {

// Establish SQL database session
Timer timer;
mysqlx::Session session{
  "mysqlx://" + modifier_password + ":" + modifier_username +
  "@localhost/" + default_database_name};

// Number of queries
uint64_t n_queries{0};

// Show an sql query
const bool show{dc("show", false)};
const std::string & show_sql(const std::string & query,
                             const std::string & placeholders = "") {
  ++n_queries;
  if (!show) return query;
  const std::string placeholders_text{maybe(
      placeholders.size(), nl + "Placeholders '?': " + placeholders)};
  std::cout << pre("Query " + str(n_queries) + ":" + nl +
                   encode(query + placeholders_text)) + nl;
  return query;
}

// Silent sql query
mysqlx::SqlResult ssql(const std::string & query) {
  return session.sql(query).execute();
}
template<class ... PHS>
mysqlx::SqlResult ssql(const std::string & query, PHS && ... placeholders) {
  return session.sql(query).bind(std::forward<PHS>(placeholders)...).execute();
}
// Possibly noisy sql query
mysqlx::SqlResult sql(const std::string & query) {
  return ssql(show_sql(query));
}
template<class ... PHS>
mysqlx::SqlResult sql(const std::string & query, PHS && ... placeholders) {
  show_sql(query, csj("'" + str(std::forward<PHS>(placeholders)) + "'" ...));
  return ssql(query, std::forward<PHS>(placeholders)...);
}

//
// Table joins
//
const StringMap key{
  {"people", "person"},
  {"campuses", "campus"},
  {"buildings", "building"},
  {"places", "place"},
  {"locations", "location"},
  {"tests", "test"},
  {"users", "user"},
  {"passwords", "password"}};
std::string table_join(const std::string & table_name,
                       const std::string & type,
                       const std::string & id = "") {
  try {
    return type + " join " + table_name + " using (" +
        (id.size() ? id : key.at(table_name)) + "_id)";
  } catch (...) {
    warning("Table name lookup failed for " + table_name);
    throw;
  }
}
std::string ljoin(const std::string & table_, const std::string & id = "") {
  return table_join(table_, "left", id);
}
std::string rjoin(const std::string & table_, const std::string & id = "") {
  return table_join(table_, "right", id);
}
std::string ijoin(const std::string & table_, const std::string & id = "") {
  return table_join(table_, "inner", id);
}

// Dates
const std::string t2020{"2020-01-01 00:00:00"};
const std::string t2100{"2100-01-01 00:00:00"};
#if 1
const std::string now{
  []() {
    time_t rawtime;
    time(&rawtime);
    struct tm * info{localtime(&rawtime)};
    static char result[20];
    strftime(result, 20, "%Y-%m-%d %H:%M:%S", info);
    return result;
  }()};
#else
const std::string now{"2021-02-28 01:01:01"};
#endif
const std::string today{now.substr(0, 10)};
std::string date_time(const std::string & col, const std::string & title_) {
  return "date_format(" + col + ", '%Y-%m-%d %H:%i:%s') as '" + title_ + "'";
}
std::string date_col(const std::string & col, const std::string & title_) {
  return "date_format(" + col + ", '%Y-%m-%d') as '" + title_ + "'";
}

//
// Queries
//
const std::string start_transaction_query{nj("start", "transaction")};
const std::string commit_transaction_query{"commit"};
const std::string rollback_transaction_query{"rollback"};
const std::string schema_query{
  nj("select", "*", "from INFORMATION_SCHEMA.TABLES",
     "where TABLE_SCHEMA = ? and TABLE_NAME = ?")};
const std::string schemas_query{
  nj("select", "*", "from INFORMATION_SCHEMA.TABLES",
     "where TABLE_SCHEMA = ?")};
std::string create_table_query(const std::string & name,
                               const std::string & definition) {
  return nj("create", "table " + name, "(", definition, ")");
}
std::string load_table_query(const std::string & name,
                             const std::string & values) {
  return nj("insert ignore", "into " + name, "values(" + values + ")");
}
const std::string optimize_query{nj("optimize", "table ?")};
const std::string show_tables_query{nj("show", "tables")};
const std::string table_columns_query{nj("show", "columns", "from ?")};
const std::string select_all_order_query{
  nj("select", "*", "from ?", "order by ?")};
const std::string foreign_query{nj("set", "session", "foreign_key_checks = ?")};
const std::string drop_table_query{nj("drop", "table ?")};

// Transaction
class Transaction {
 public:
  explicit Transaction(const std::string & name_ = "") : name{name_} {
    sql(start_transaction_query);
  }
  void commit() {
    sql(commit_transaction_query);
    committed = true;
  }
  ~Transaction() {
    if (!committed) {
      sql(rollback_transaction_query);
      warning("Transaction " + name + " rolled back");
    }
  }

 private:
  std::string name{};
  bool committed{false};
};

// Check if sql command affected just one row
bool one(const mysqlx::SqlResult & result, const std::string & message) {
  if (result.getAffectedItemsCount() == 1) {
    success("Success: " + message);
    return true;
  } else {
    warning("Failure: " + message);
    return false;
  }
}

// Foreign key def
std::string foreign(const std::string & table_,
                    const std::string & name = "") {
  try {
    const std::string col{key.at(table_) + "_id"};
    const std::string fk{name.size() ? name + "_id" : col};
    return "foreign key (" + fk + ") references " + table_ + "(" + col + ")";
  } catch (...) {
    warning("Key for table not found: " + table_);
    throw;
  }
}

// Load passed values into database
template <class Fun>
void load(const std::string &, const Fun & function = []() {}) {
  function();
}
// Load passed values into database (never pass untrusted values)
void load(const std::string & name, const Strings & values) {
  for (const std::string & value : values)
    sql(load_table_query(name, value));
}

// Fill placeholders (do not use with user data)
std::string fph(const std::string & text) { return text; }
template <class String, class ... PHS>
std::string fph(std::string text, String && first, PHS && ... phs) {
  return fph(replace_substring(text, " ?",
                               std::string(" ") + std::forward<String>(first)),
              std::forward<PHS>(phs) ...);
}

// Populate a table with canned info
Strings table_names;
template <class Load>
void table(const std::string & name, const std::string & definition,
           const Load & values) {
  try {
    // Create table if nonexistent
    if ((reset || reload) &&
        sql(schema_query, default_database_name, name).count() == 0) {
      sql(create_table_query(name, definition));
      load(name, values);
      sql(fph(optimize_query, name));
    }

    // Return name
    table_names.push_back(name);
  } catch (const mysqlx::Error & e) {
    warning("Problem getting table " + name);
    throw;
  }
}

// Reset entire database!
void reset_database() {
  sql(foreign_query, 0);
  for (const mysqlx::Row row : sql(show_tables_query)) {
    const std::string table_name{std::string(row[0])};
    if (table_name != "backups")
      sql(fph(drop_table_query, table_name));
  }
  sql(foreign_query, 1);
  success("Database was cleared and repopulated with basic info");
}
const mysqlx::Value null{nullptr};

// Show any RowResult or SqlResult object as a table
uint64_t tid{0};
std::string table_id;
const bool csv{dc("csv", false)};
const bool scroll{dc("scroll", true)};
template <class ROWS>
uint64_t table(const std::string & title_, ROWS && rows,
               const std::string & extra = "") {
  std::ostringstream out;
  const uint64_t n_columns{rows.getColumnCount()};
  out << h2(title_ + " " + nl +
            span("(" + str(rows.count()) + "&nbsp;entr" +
                 (rows.count() == 1 ? "y" : "ies") + ")" +
                 (extra.size() ? "</span> " + nl + "<span>" + extra : ""),
                 {{"class", "entries"}}),
                 anchor("h" + (table_id = str(++tid)))) << nl;
  const uint64_t n_rows{rows.count()};
  if (!n_rows) {
    std::cout << out.str();
    warning("There was no data to show");
    return n_rows;
  }
  if (csv) {
    out << "<pre>" << nl;
    for (uint64_t c{0}; c != n_columns; ++c)
      out << maybe(c, ",") <<
          replace_all(std::string(rows.getColumn(c).getColumnLabel()), br, " ");
    out << nl;
    for (const mysqlx::Row row : rows) {
      for (uint64_t c{0}; c != n_columns; ++c) out << maybe(c, ",") << row[c];
      out << nl;
    }
    out << "</pre>" << nl;
  } else {
    out << table_tag({{"class", "sql" + maybe(scroll, " scroll")},
                      {"id", "t" + table_id}}, false) << nl << "<tr>" << nl;
    for (uint64_t c{0}; c != n_columns; ++c)
      out << th(rows.getColumn(c).getColumnLabel()) << nl;
    out << th("Row") << nl << "</tr>" << nl;
    for (uint64_t r{0}; r != n_rows ; ++r) {
      const mysqlx::Row row{rows.fetchOne()};
      out << "<tr>" << nl;
      for (uint64_t c{0}; c != n_columns; ++c)
        out << "<td>" << row[c] << "</td>" << nl;
      out << "<td>" << r + 1 << "</td>" << nl << "</tr>" << nl;
    }
    out << "</table>" << nl;
  }
  std::cout << out.str();
  return n_rows;
}

template<class Result>
void spanify(const std::string & text, Result && result) {
  table(replace_one(text, ":", ": " + nl + "<span>") + "</span>", result);
}
void spanify(const std::string & text, const std::string & query) {
  mysqlx::SqlResult result{sql(query)};
  if (result.count()) spanify(text, result);
}
template<>
void spanify(const std::string & text, std::string && query) {
  spanify(text, static_cast<const std::string &>(query));  // Why cast needed?
}

const Strings sensitive_tables{
  "tests", "backups", "queries", "users", "passwords"};

// Show any table
void table(const std::string & title_, mysqlx::Table & table_,
           const bool permission, const bool show_schema = true) {
  const std::string name{table_.getName()};
  if (show_schema)
    table("Schema for " + name, sql(fph(table_columns_query, name)));
  if (permission == false) {
    for (const std::string & name_ : sensitive_tables) {
      if (name == name_) {
        warning("The " + name_ +
                " table requires additional privileges to view");
        return;
      }
    }
  }
  const std::string order(sql(fph(table_columns_query, name)).fetchOne()[0]);
  mysqlx::SqlResult result{sql(fph(select_all_order_query, name, order))};
  table(title_, result);
}

struct ColumnDef {
  ColumnDef(const std::string & name_, const std::string & def_) :
      name{name_}, def{def_} {}
  const std::string name;
  const std::string def;
  bool primary{false};
  bool key{false};
  std::unique_ptr<StringPair> foreign{nullptr};
};
struct TableDef {
  using Columns = std::vector<ColumnDef>;
  explicit TableDef(const std::string & name_) : name{name_} {}
  const std::string name;
  uint64_t size() const { return columns.size(); }
  ColumnDef & operator[](const std::string & column) {
    try {
      return columns[column_lookup.at(column)];
    } catch (...) {
      throw Error("Bad column lookup for") << column << "in" << name;
    }
  }
  const ColumnDef & operator[](const uint64_t column) const {
    return columns[column];
  }
  const ColumnDef & operator[](const std::string & column) const {
    try {
      return columns[column_lookup.at(column)];
    } catch (...) {
      throw Error("Bad column lookup for") << column << "in" << name;
    }
  }
  void add(const std::string & column_name, const std::string & column_def) {
    if (column_lookup.emplace(column_name, columns.size()).second == false)
      warning("Duplicate column " + column_name + " for table " + name);
    columns.emplace_back(column_name, column_def);
  }
  Columns::iterator begin() { return columns.begin(); }
  Columns::iterator end() { return columns.end(); }
  Columns::const_iterator begin() const { return columns.begin(); }
  Columns::const_iterator end() const { return columns.end(); }

 private:
  Columns columns{};
  std::map<std::string, uint64_t> column_lookup{};
};
using TableDefs = std::map<std::string, TableDef>;
std::string table_create_query(const std::string & name) {
  return fph("show create table ?", name);
}
std::string table_regex_string{
  R"xxx(CREATE TABLE `(.+)`\s+\(\s+([\s\S]+)\s+\)\s+(ENGINE.+)\n?)xxx"};
std::regex line_regex{R"xxx((.*?)\s?`(\w+)`\s?(.+))xxx"};
std::regex table_regex{table_regex_string};
std::string database_diagram(const Strings & table_order =
                             {"queries", "locations", "people", "bosses",
                              "campuses", "buildings", "places", "tests",
                              "passwords", "users", "calendar",  "layouts",
                              "", "", "", "backups"}) {
  std::ostringstream out;
  const bool verbose{false};
  TableDefs tables;
  for (const mysqlx::Row row : sql(show_tables_query)) {
    const std::string table_name{std::string(row[0])};
    TableDef table{table_name};
    const std::string table_command{
      std::string(sql(table_create_query(table_name)).fetchOne()[1])};
    std::smatch table_matches;
    if (regex_match(table_command, table_matches, table_regex) == false)
      warning("Bad regex search on table " + table_name);
    if (verbose) out << h3(table_name) << nl;
    if (verbose) out << pre(table_command) << nl;
    for (const std::string & line : split(table_matches[2].str(), ",\n  ")) {
      std::string name;
      std::string def;
      std::istringstream line_stream{line};
      switch (line[0]) {
        case '`':
          {
            line_stream.get();
            getline(line_stream, name, '`');
            line_stream.get();
            getline(line_stream, def);
            if (verbose) out << p("Column " + name + " = " + def) << nl;
            table.add(name, def);
          }
          break;
        case 'P':
          {
            getline(line_stream, def, '(');
            if (def != "PRIMARY KEY ") warning("Bad primary key def " + def);
            getline(line_stream, def, ')');
            if (verbose) out << p("Primary key " + def) << nl;
            const Strings keys{split(def, ",")};
            for (std::string key_name : keys) {
              key_name = key_name.substr(1, key_name.size() - 2);
              table[key_name].primary = true;
              table[key_name].key = true;
            }
          }
          break;
        case 'K':
        case 'U':
          {
            getline(line_stream, def, '`');
            if (def.find("KEY ") == std::string::npos)
              warning("Bad key def " + def + " for " + line);
            std::string key_type{def.find("UNIQUE") != std::string::npos ?
              "Unique Key" : "Key"};
            if (key_type == "Key" && def != "KEY ")
              warning("Uncovered key type " + def);
            getline(line_stream, name, '`');
            getline(line_stream, def, '(');
            getline(line_stream, def, ')');
            if (verbose) out << p(key_type + " " + name + " = " + def) << nl;
            const Strings keys{split(def, ",")};
            for (std::string key_name : keys) {
              key_name = key_name.substr(1, key_name.size() - 2);
              table[key_name].key = true;
            }
          }
          break;
        case 'C':
          {
            getline(line_stream, def, '`');
            if (def != "CONSTRAINT ") warning("Bad constraint for " + line);
            getline(line_stream, name, '`');
            getline(line_stream, def, '`');
            if (def != " FOREIGN KEY (") warning("Bad foreign for " + line);
            std::string local_column;
            std::string foreign_column;
            std::string foreign_table;
            getline(line_stream, local_column, '`');
            getline(line_stream, def, '`');
            if (def != ") REFERENCES ") warning("Bad references for " + line);
            getline(line_stream, foreign_table, '`');
            getline(line_stream, def, '`');
            getline(line_stream, foreign_column, '`');
            if (verbose) out << p("Constraint " + name + " : " + local_column +
                                  " references " + foreign_column + " in " +
                                  foreign_table) << nl;
            table[local_column].foreign =
                std::make_unique<StringPair>(foreign_table, foreign_column);
          }
          break;
        default:
          warning("Unknown line type X" + line);
      }
      if (!line_stream) warning("Bad parse for line " + line);
    }
    tables.emplace(table_name, std::move(table));
  }
  // Check constraints
  uint64_t n_constraints{0};
  for (const auto & table_info : tables) {
    const TableDef & table{table_info.second};
    for (const ColumnDef & column : table) {
      const std::unique_ptr<StringPair> & foreign{column.foreign};
      if (foreign) {
        ++n_constraints;
        const std::string table_name{foreign->first};
        const std::string column_name{foreign->second};
        auto found = tables.find(table_name);
        if (found == tables.end()) {
          warning("Missing foreign table " + table_name + " for column " +
                  column.name + " in " + table.name);
        } else {
          const TableDef & foreign_table{found->second};
          foreign_table[column_name];
        }
      }
    }
  }
  // Layout diagram
  const uint64_t n_columns{4};
  using Coords = std::pair<uint64_t, uint64_t>;
  std::map<std::string, Coords> lookup;
  uint64_t n{0};
  Strings names;
  for (const std::string & name : table_order) {
    // const TableDef & table{tables[name]};
    if (name.size()) {
      lookup[name] = Coords{n / n_columns, n % n_columns};
      names.push_back(name);
    }
    ++n;
  }
  auto distance = [](const Coords & l, const Coords & r) {
    const double misorder_penalty{0.5};
    return (l.first > r.first ? l.first - r.first :
            r.first - l.first + misorder_penalty) +
        1.5 * (l.second > r.second ? l.second - r.second :
               r.second - l.second + misorder_penalty);
  };
  using Total = std::pair<double, std::string>;
  auto get_length = [&tables, &lookup, &distance]() {
    double total_distance{0};
    double worst_distance{0};
    std::string worst_node{""};
    for (const auto & table_info : tables) {
      const TableDef & table{table_info.second};
      const std::string table_name{table.name};
      const Coords coords{lookup[table_name]};
      for (const ColumnDef & column : table) {
        const std::unique_ptr<StringPair> & foreign{column.foreign};
        if (foreign) {
          const std::string foreign_table_name{foreign->first};
          const Coords foreign_coords{lookup[foreign_table_name]};
          const double link_distance{distance(coords, foreign_coords)};
          total_distance += link_distance;
          if (link_distance > worst_distance) {
            worst_distance = link_distance;
            worst_node = table_name;
          }
        }
      }
    }
    return Total{total_distance, worst_node};
  };
  if (0) {
    const Total first_try{get_length()};
    double best_length{first_try.first};
    double last_length{best_length * 2};
    while (true) {
      if (verbose) out << p(std::to_string(best_length) + " " +
                            std::to_string(n_constraints)) << nl;
      last_length = best_length;
      StringPair best_swap{"", ""};
      for (const std::string & name1 : names) {
        for (const std::string & name2 : names) {
          if (name1 == name2) continue;
          std::swap(lookup[name1], lookup[name2]);
          const Total this_try{get_length()};
          if (this_try.first < best_length) {
            best_length = this_try.first;
            best_swap = {name1, name2};
          }
          std::swap(lookup[name1], lookup[name2]);
        }
      }
      if (best_length < last_length) {
        std::swap(lookup[best_swap.first], lookup[best_swap.second]);
      } else {
        break;
      }
    }
  }

  // Make svg
  std::map<Coords, std::string> reverse;

  uint64_t last_y{100};
  for (const auto & info : lookup) reverse.emplace(info.second, info.first);
  const double font_size{14};
  const double cell_border{1};
  const double cell_padding{2};
  const double table_x_spacing{30};
  const double table_y_spacing{20};
  const double table_x_margin{table_y_spacing / 2};
  const double table_y_margin{10};
  const double cell_text_width{font_size * 16};
  const double cell_height{font_size + 2 * (cell_padding + cell_border)};
  const double cell_width{cell_text_width + 2 * (cell_padding + cell_border)};
  const double table_x_offset{cell_width + table_x_spacing};
  const double svg_width{2 * table_x_margin + n_columns * cell_width +
    (n_columns - 1) * table_x_spacing};
  std::vector<double> heights(n_columns, table_y_margin);
  std::ostringstream svg;
  using dCoords = std::pair<double, double>;
  std::map<std::string, dCoords> col_coords;
  std::vector<StringPair> links;
  mysqlx::Schema db{session.getSchema(default_database_name)};
  for (const auto & info : reverse) {
    const std::string table_name{info.second};
    const uint64_t n_rows{db.getTable(table_name).count()};
    const TableDef & table{tables.at(table_name)};
    const Coords coords{info.first};
    const uint64_t row{coords.first};
    const uint64_t col{coords.second};
    const double x{table_x_margin + col * table_x_offset};
    const double y{heights[col]};
    const uint64_t n_lines{table.size() + 1};
    const double table_height{n_lines * cell_height};
    svg << "<g>" << nl;
    svg << "<rect x=\"" << x << "\" y=\"" << y << "\""
        << " width=\"" << cell_width << "\""
        << " height=\"" << cell_height << "\" />\n"
        << "<a xlink:href=\"?view=Tables&name=" << table_name << "\">"
        << "<text x=\"" << x + cell_width / 2 << "\""
        << " y=\"" << y + cell_height / 2 << "\">"
        << table_name << "<title>" << n_rows << " rows</title></text></a>\n";
    for (uint64_t c{0}; c != table.size(); ++c) {
      const ColumnDef & column{table[c]};
      const std::string column_name{column.name};
      const double ly{y + (c + 1) * cell_height};
      const std::string name{table_name + "." + column_name};
      col_coords.emplace(name, dCoords{ly + 0.5 * cell_height, x});
      if (column.foreign) {
        links.emplace_back(
            name, column.foreign->first + "." + column.foreign->second);
      }
      const std::string class_string{
        column.key ? " class=\"key\"" : ""};
      svg << "<rect x=\"" << x << "\" y=\"" << ly << "\""
          << " width=\"" << cell_width << "\""
          << " height=\"" << cell_height << "\" />\n"
          << "<text" << class_string
          << " x=\"" << x + cell_padding + cell_border << "\""
          << " y=\"" << ly + cell_height / 2 << "\">"
          << column_name << (column.primary ? "*" : "") << "</text>\n"
          << "<text class=\"right\""
          << " x=\"" << x + cell_padding + cell_border + cell_text_width << "\""
          << " y=\"" << ly + cell_height / 2 << "\">"
          << first_word(column.def)
          << "<title>" << column.def << "</title></text>\n";
    }
    svg << "</g>" << nl;
    if (row != last_y) last_y = row;
    heights[col] += table_height + table_y_spacing;
  }
  // Make links
  for (const StringPair & link : links) {
    const std::string from{link.first};
    const std::string to{link.second};
    const dCoords from_coords{col_coords[from]};
    const dCoords to_coords{col_coords[to]};
    const bool right_from{from_coords.second < to_coords.second};
    const bool right_to{to_coords.second < from_coords.second};
    const double from_x{from_coords.second + (right_from ? cell_width : 0.0)};
    const double to_x{to_coords.second + (right_to ? cell_width : 0.0)};
    const double from_y{from_coords.first};
    const double to_y{to_coords.first};
    // const double mid_x{(from_x + to_x) / 2};
    const double mid_y{(from_y + to_y) / 2};
    const std::string temp{" title=\"" + from + "->" + to + "\""};
    if (from_coords.second < to_coords.second ||
        from_coords.second > to_coords.second) {
      svg << "<path" << temp << " d=\"M" << from_x << " " << from_y
          << "L" << to_x << " " << to_y << "\"/>" << nl;
    } else {
      const double d{table_x_spacing / 2};
      svg << "<path" << temp << " d=\"M" << from_x << " " << from_y
          << "Q" << from_x - d << " " << from_y << " "
          << from_x - d << " " << mid_y << " "
          << "Q" << from_x - d << " " << to_y << " "
          << to_x << " " << to_y << "\"/>" << nl;
    }
  }
  const double svg_height{*max_element(heights.begin(), heights.end()) -
    table_y_spacing + table_y_margin};
  const double padding{100.0 * svg_height / svg_width};
  std::ostringstream svg_elem;
  svg_elem << nj(
      "<div class=\"db\" style=\"width:100%; padding-top:" +
      str(padding) + "%;\">",
      "<svg viewBox=\"0 0 " + str(svg_width) + " " + str(svg_height) + "\">",
      svg.str() + "</svg>", "</div>") << nl;
  if (verbose) {
    svg_elem << pre(encode(svg_elem.str())) << nl;
    svg_elem << out.str();
  }
  return svg_elem.str();
}

}  // namespace paa

#endif  // UTILITY_SQL_H_
