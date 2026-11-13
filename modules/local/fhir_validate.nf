include { asFlag } from '../../subworkflows/local/utils_params'

process FHIR_VALIDATE {
    label 'process_medium'
    conda "conda-forge::openjdk=21"
    container "eclipse-temurin:21-jre"

    input:
    path fhir_files, stageAs: 'input/*'
    path validator_jar
    path codesystems

    output:
    tuple path("input/*.json", includeInputs: true), path("validation.json"), emit: result
    path "validator.log",                                                    emit: log
    path "versions.yml",                                                     emit: versions

    script:
    def mem = (task.memory.toMega() * 0.8).intValue()
    """
    fhir_home="\${HOME:-}"
    if [ -z "\$fhir_home" ] || [ ! -w "\$fhir_home" ]; then
        fhir_home="\$PWD/fhir_home"
        mkdir -p "\$fhir_home"
    fi

    java -Xmx${mem}M -Dfile.encoding=UTF-8 -Duser.home="\$fhir_home" -jar ${validator_jar} \\
        -version 4.0.1 \\
        -ig ${params.fhir_ig} \\
        -ig ${codesystems} \\
        -tx ${params.fhir_tx_server} \\
        -output validation.json \\
        input/*.json > validator.log 2>&1 || true

    if [ ! -s validation.json ]; then
        tail -50 validator.log >&2
        echo "ERROR: the FHIR validator did not produce output" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fhir_validator: \$(sed -n 's/.*FHIR Validation tool Version \\([0-9.]*\\).*/\\1/p' validator.log | head -1)
        implementation_guide: ${params.fhir_ig}
    END_VERSIONS
    """

    stub:
    """
    echo '{"resourceType": "OperationOutcome", "issue": []}' > validation.json
    touch validator.log
    printf '"%s":\\n    fhir_validator: stub\\n' "${task.process}" > versions.yml
    """
}

process FHIR_VALIDATION_SUMMARY {
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple path(fhir_files, stageAs: 'input/*'), path(validation_json)

    output:
    path "validation_summary.tsv", emit: summary
    path "passed/*",               emit: passed, optional: true
    path "versions.yml",           emit: versions

    script:
    def configured = asFlag(params.fail_on_validation_error)
    def fail_on_error = configured != null ? configured : params.fhir_profile == 'generic'
    def args = []
    if (fail_on_error) args << '--fail_on_error'
    if (params.fhir_tx_server == 'n/a') args << '--offline'
    """
    summarize_validation.py --validation ${validation_json} --files input/*.json ${args.join(' ')}

    mkdir -p passed
    while read -r name; do
        [ -n "\$name" ] && cp "input/\$name" passed/
    done < passed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'file\\terrors\\twarnings\\tinformation\\n' > validation_summary.tsv
    for f in input/*.json; do printf '%s\\t0\\t0\\t0\\n' "\$(basename "\$f")" >> validation_summary.tsv; done
    mkdir passed
    cp input/*.json passed/ 2>/dev/null || true
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
