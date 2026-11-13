include { validate } from 'plugin/nf-schema'

def schemaTypes(String schema_path) {
    def schema = new groovy.json.JsonSlurper().parse(file(schema_path))
    def types = [:]
    schema['$defs'].each { _group, defn -> defn.properties.each { name, spec -> types[name] = spec.type } }
    return types
}

def typedValue(type, value) {
    if (type == null || value == true || value == false) {
        return value
    }
    def text = value.toString().trim()
    if (type == 'integer' && text ==~ /-?\d{1,18}/) {
        def number = text.toLong()
        return number >= Integer.MIN_VALUE && number <= Integer.MAX_VALUE ? number.toInteger() : number
    }
    if (type == 'number' && text.isNumber()) {
        return text.toBigDecimal()
    }
    if (type == 'boolean' && text.toLowerCase() in ['true', 'false']) {
        return text.toLowerCase() == 'true'
    }
    if (type == 'string') {
        return value.toString()
    }
    return value
}

def validateTypedParams() {
    def schema_path = "${projectDir}/nextflow_schema.json".toString()
    def types = schemaTypes(schema_path)
    def unknown = params.keySet().findAll { name -> !types.containsKey(name) }
    if (unknown) {
        log.warn "Unknown parameters, ignored: " + unknown.collect { name -> "--${name}" }.join(', ')
    }
    def typed = [:]
    params.findAll { _name, value -> value != null }.each { name, value -> typed[name] = typedValue(types[name], value) }
    try {
        validate(typed, schema_path)
    } catch (Exception e) {
        def failure = e.cause ?: e
        def lines = (failure.message ?: failure.toString()).readLines()
        def problems = lines.findAll { line -> line.trim().startsWith('*') || line.trim().startsWith('/') }
        def report = problems ? problems.collect { line -> '  ' + line.trim().replaceFirst(/^\* */, '').replaceFirst(/^\//, '--') } : lines
        error("Invalid parameters:\n" + report.join('\n'))
    }
}

def asFlag(value) {
    if (value == null) {
        return null
    }
    def text = value.toString().trim().toLowerCase()
    if (text in ['true', 'false']) {
        return text == 'true'
    }
    error("Expected true or false, got '${value}'")
}

def asInt(value) {
    if (value == null) {
        return null
    }
    def text = value.toString().trim()
    if (text ==~ /-?\d{1,9}/) {
        return text.toInteger()
    }
    error("Expected an integer, got '${value}'")
}
