
package us.kbase.kbvariationimporter;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: import_snp_params</p>
 * <pre>
 * Insert your typespec information here.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "staging_file_subdir_path",
    "will_perform_gwas"
})
public class ImportSnpParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("staging_file_subdir_path")
    private String stagingFileSubdirPath;
    @JsonProperty("will_perform_gwas")
    private Long willPerformGwas;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public ImportSnpParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("staging_file_subdir_path")
    public String getStagingFileSubdirPath() {
        return stagingFileSubdirPath;
    }

    @JsonProperty("staging_file_subdir_path")
    public void setStagingFileSubdirPath(String stagingFileSubdirPath) {
        this.stagingFileSubdirPath = stagingFileSubdirPath;
    }

    public ImportSnpParams withStagingFileSubdirPath(String stagingFileSubdirPath) {
        this.stagingFileSubdirPath = stagingFileSubdirPath;
        return this;
    }

    @JsonProperty("will_perform_gwas")
    public Long getWillPerformGwas() {
        return willPerformGwas;
    }

    @JsonProperty("will_perform_gwas")
    public void setWillPerformGwas(Long willPerformGwas) {
        this.willPerformGwas = willPerformGwas;
    }

    public ImportSnpParams withWillPerformGwas(Long willPerformGwas) {
        this.willPerformGwas = willPerformGwas;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((("ImportSnpParams"+" [workspaceName=")+ workspaceName)+", stagingFileSubdirPath=")+ stagingFileSubdirPath)+", willPerformGwas=")+ willPerformGwas)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
