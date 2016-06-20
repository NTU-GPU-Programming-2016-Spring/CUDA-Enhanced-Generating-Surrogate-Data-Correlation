// Native libraries below.
var fs = require('graceful-fs');
// 3rd-party libraries below.
var _  = require('underscore');

var csvPath    = __dirname + '/../data/original/';
var csvNewPath = __dirname + '/../data/uncompressed/';

fs.readdir(csvPath, function (err, files) {
	// Each file.
	files.map(function (file) {
		// Read the file.
		fs.readFile(csvPath + file, 'utf8', function (err, data) {
			// Split by new line.
			var rowList = data.split('\n');
			rowList = _.filter(rowList, function (str) { return str != '' });
			// Write into CSV each row.
			rowList.map(function (row, index) {
				var filename = file.match(/^(.+)\.csv$/)[1] + '-' + (index + 1) + '.csv';
				fs.writeFile(csvNewPath + filename, row, function (err) { });
			});
		});
	});
});