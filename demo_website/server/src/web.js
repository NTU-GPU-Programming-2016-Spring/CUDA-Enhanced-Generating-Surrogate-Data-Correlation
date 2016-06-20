var exports = module.exports = {};

// Native libraries below.
var fs     = require('fs');
// 3rd-party libraries below.
var _      = require('underscore');
var async  = require('async');
var jade   = require('jade');

// Static variables.
var subjects = {
	'AD1':  [1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1],
	'SUP1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
	'AD2':  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
	'SUP2': [0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0]
}

exports.core = {
	index: function * () {
		var html  = jade.renderFile('./views/index.jade', {});
		this.body = html;
	}
};

exports.chart = {
	update: function * () {
		var body = this.request.body;
		if (! body.position) return;
		var json = {
			'cc': yield parseCCData(body.view, body.brain, body.position),
			'subjectListHTML' : jade.renderFile('./views/subject.jade', { 'subjectList': subjects[body.view] })
		};
		this.type = 'application/json';
		this.body = json;
	}
};

function parseCCData(view, brain, pos) {
	var filename = __dirname + '/../data/uncompressed/2000-' + view + '-' + (brain == 'left' ? 'lh' : 'rh') + '-' + pos + '.csv';

	return new Promise(function (resolve, reject) {
		fs.readFile(filename, 'utf8', function (err, data) {
			if (err) return resolve([]);
			var arrData = data.split(',').map(function (str) { return parseFloat(parseFloat(str).toFixed(4)); });
			return resolve(arrData);
		});
	});
}